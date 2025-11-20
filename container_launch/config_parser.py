#!/usr/bin/env python3
"""
plan.py - config loader / merger / validator / reporter for RUFUS.

Usage (examples):
  # basic plan (use default_profile if present)
  ./plan.py --config /work/run-config.yaml --workdir /work

  # explicitly choose a profile
  ./plan.py --config /work/run-config.yaml --profile whole_genome_default --workdir /work

  # validate only (exit non-zero on problems)
  ./plan.py --config /work/run-config.yaml --validate --workdir /work

  # describe a profile (what fields it introduces / requires)
  ./plan.py --config /work/run-config.yaml --describe-profile windowed_slurm

  # generate a minimal template for a profile
  ./plan.py --config /work/run-config.yaml --generate-template windowed_slurm > my-template.yaml

  # dump the fully resolved config to a file
  ./plan.py --config /work/run-config.yaml --dump-effective resolved.yaml --workdir /work
"""

import argparse
import copy
import os
import sys
import textwrap
import yaml
import json

# ---------- utilities ----------
def deep_merge(a, b):
    """
    Recursively merge dict b into dict a and return a new dict.
    Values in b override values in a. Lists are replaced, not merged.
    """
    if a is None:
        a = {}
    out = copy.deepcopy(a)
    for k, v in (b or {}).items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out

def load_yaml(path):
    with open(path, 'r') as fh:
        return yaml.safe_load(fh) or {}

def dump_yaml(obj, fh):
    yaml.safe_dump(obj, fh, sort_keys=False)

def fmt_bool(b):
    return 'true' if b else 'false'

def print_header(msg):
    print("\n" + "="*len(msg))
    print(msg)
    print("="*len(msg))

# ---------- schema helpers (minimal doc) ----------
# Minimal, human-oriented schema info used by describe_profile/generate_template
SCHEMA_HINTS = {
    "inputs.subject_file":       ("path", "Required: subject sample (CRAM/BAM)"),
    "inputs.paired_control_files": ("array", "Optional: list of control CRAM/BAM, used if method == paired_files"),
    "inputs.reference_fasta":    ("path", "Required: reference FASTA"),
    "controls.method":           ("enum", "Control selection method: 'prebuilt' or 'paired_files'"),
    "prebuilt.control.control_hash_version": ("string", "Default control hash version to fetch"),
    "prebuilt.control.control_hash_local_dir": ("path", "Optional path to local control hashes"),
    "prebuilt.kg1.kg1_hash_version": ("string", "Default 1KG hash version to fetch"),
    "prebuilt.kg1.kg1_hash_local_dir": ("path", "Optional path to local 1KG hashes"),
    "window.window_size":        ("int", "Window size in bp (region mode) - supports only 1000 for now"),
    "window.region_file":        ("path", "Region file required for windowed runs"),
    "tools.container_type":      ("enum", "docker or singularity"),
    "tools.resource_mgmt":       ("enum", "slurm or parallel"),
    "container.singularity_image":("path", "Path to SIF image (if using singularity)"),
    "container.docker_image":    ("string", "Docker image name (if using docker)"),
    "resources.slurm.account":   ("string", "Slurm account (required for slurm jobs)"),
    "resources.slurm.partition": ("string", "Slurm partition"),
    "resources.slurm.time":      ("string", "Slurm walltime like DD-HH:MM:SS"),
    "resources.slurm.email":     ("string", "Email for Slurm notifications (optional)"),
}

# ---------- resolver / validator ----------
class Planner:
    def __init__(self, cfg, cli_profile=None, workdir=None, cli_overrides=None):
        self.raw_cfg = cfg or {}
        self.cli_profile = cli_profile
        self.workdir = os.path.abspath(workdir or os.getcwd())
        self.cli_overrides = cli_overrides or {}
        self.selected_profile_name = None
        self.profile_cfg = {}
        self.effective = {}
        self.warnings = []
        self.errors = []

    def select_profile(self):
        profiles = self.raw_cfg.get('profiles', {}) or {}
        # precedence: CLI profile > raw_cfg.default_profile > None
        p = self.cli_profile or self.raw_cfg.get('default_profile')
        if p:
            if p not in profiles:
                self.errors.append(f"Profile '{p}' not found in config (available: {list(profiles.keys())})")
                return
            self.selected_profile_name = p
            self.profile_cfg = profiles.get(p) or {}
        else:
            # no profile defined - empty profile
            self.selected_profile_name = None
            self.profile_cfg = {}

    def merge(self):
        base = copy.deepcopy(self.raw_cfg)
        # remove profiles block from base so it does not appear in effective
        base.pop('profiles', None)
        base.pop('default_profile', None)
        merged = deep_merge(base, self.profile_cfg)
        # apply CLI overrides (explicit simple keys allowed)
        if self.cli_overrides:
            merged = deep_merge(merged, self.cli_overrides)
        self.effective = merged

    def detect_overrides(self):
        # compare top-level raw_cfg (without profiles) vs profile_cfg for overlapping keys
        top = copy.deepcopy(self.raw_cfg)
        top.pop('profiles', None)
        top.pop('default_profile', None)
        # walk keys from profile and see if they exist in top at same path and differ
        def walk(d, prefix=[]):
            for k, v in (d or {}).items():
                path = prefix + [k]
                yield path, v
                if isinstance(v, dict):
                    yield from walk(v, path)
        for path, pval in walk(self.profile_cfg):
            # try to get value at same path in top
            cur = top
            found = True
            for part in path:
                if isinstance(cur, dict) and part in cur:
                    cur = cur[part]
                else:
                    found = False
                    break
            if found:
                if cur != pval:
                    self.warnings.append(f"PROFILE OVERRIDE: profile '{self.selected_profile_name}' overrides top-level {'.'.join(path)}: {cur!r} -> {pval!r}")

    def resolve_prebuilt(self):
        """
        Resolve whether control and kg1 prebuilt hashes are local or will be fetched.
        Add resolved keys under self.effective['_resolved']['prebuilt'][...]
        """
        prebuilt = self.effective.get('prebuilt', {}) or {}
        resolved = {'control': {}, 'kg1': {}}
        for key, cfg in (('control', prebuilt.get('control', {}) or {}), ('kg1', prebuilt.get('kg1', {}) or {})):
            local_dir = cfg.get(f'{ "control_hash_local_dir" if key=="control" else "kg1_hash_local_dir" }')
            # keys differ in name; normalize:
            if key == 'control':
                version = cfg.get('control_hash_version')
                local = cfg.get('control_hash_local_dir') or ''
            else:
                version = cfg.get('kg1_hash_version')
                local = cfg.get('kg1_hash_local_dir') or ''
            if local:
                resolved[key]['use_local'] = True
                resolved[key]['local_dir'] = local
                resolved[key]['fetch_version'] = None
            else:
                resolved[key]['use_local'] = False
                resolved[key]['local_dir'] = None
                # set a recommended download location under workdir
                if version:
                    out = os.path.join(self.workdir, 'prebuilt', key, version)
                else:
                    out = os.path.join(self.workdir, 'prebuilt', key, 'unknown_version')
                resolved[key]['fetch_version'] = version
                resolved[key]['recommended_local_path'] = out
        self.effective.setdefault('_resolved', {})['prebuilt'] = resolved

    # --- helper to detect common placeholders ---
    # --- helper to detect common placeholders (if not already present) ---
def is_placeholder(s):
    if s is None:
        return False
    if not isinstance(s, str):
        return False
    low = s.lower()
    placeholders = (
        "your_", "your-", "youraccount", "your_account", "your_partition",
        "path_to/", "path_to\\", "example", "replace", "your_email", "your_unid",
        "your_account_here", "your_partition_here"
    )
    return any(tok in low for tok in placeholders)

# --- Replace Planner.validate with this function ---
def validate(self):
    """
    Strong validation rules. Appends fatal errors to self.errors and
    non-fatal issues to self.warnings.
    Enforces:
      - required top-level inputs exist and point to real files
      - region mode requires window_size or a region file that exists
      - controls.method in ('prebuilt','paired_files','both')
      - paired files present and exist in required modes (fatal)
      - prebuilt config includes version or local_dir (fatal)
      - KG1 config required unless no_kg1_removal is true (fatal)
      - slurm/container checks as before
    """
    cfg = self.effective

    # 1) Required top-level inputs presence + existence
    inputs = (cfg.get('inputs') or {})
    subj = inputs.get('subject_file')
    ref = inputs.get('reference_fasta')

    if not subj:
        self.errors.append("Missing required field: inputs.subject_file")
    else:
        if is_placeholder(subj):
            self.warnings.append("inputs.subject_file looks like a placeholder; replace with a real path.")
        elif not os.path.exists(subj):
            self.errors.append(f"Subject file not found: {subj}")

    if not ref:
        self.errors.append("Missing required field: inputs.reference_fasta")
    else:
        if is_placeholder(ref):
            self.warnings.append("inputs.reference_fasta looks like a placeholder; replace with a real path.")
        elif not os.path.exists(ref):
            self.errors.append(f"Reference FASTA not found: {ref}")

    # 2) threads_per_job sanity
    tpj = cfg.get('threads_per_job')
    try:
        if tpj is None:
            self.warnings.append("threads_per_job not set; using pipeline default.")
        else:
            if int(tpj) <= 0:
                self.errors.append("threads_per_job must be a positive integer")
    except Exception:
        self.errors.append("threads_per_job must be an integer")

    # 3) controls.method must exist and be valid
    method = (cfg.get('controls') or {}).get('method')
    if method is None:
        self.errors.append("controls.method must be set in the effective config (use 'prebuilt', 'paired_files', or 'both').")
        return
    if method not in ('prebuilt', 'paired_files', 'both'):
        self.errors.append(f"controls.method has invalid value: {method!r}. Allowed: 'prebuilt', 'paired_files', 'both'")

    # helper: collect paired files list (profile overrides controls.paired_files; fallback to inputs)
    pf = None
    controls_block = cfg.get('controls') or {}
    if isinstance(controls_block.get('paired_files'), list) and controls_block.get('paired_files'):
        pf = controls_block.get('paired_files')
    else:
        pf = (cfg.get('inputs') or {}).get('paired_control_files') or []

    # 4) Enforce paired-files depending on method (per your requested policy)
    # You requested that for method == 'paired_files' OR method == 'prebuilt' OR method == 'both',
    # missing or non-existent paired control files should be FATAL.
    if method in ('paired_files', 'prebuilt', 'both'):
        # require non-empty list
        if not isinstance(pf, list) or len(pf) == 0:
            self.errors.append(f"controls.method == '{method}' requires inputs.paired_control_files or controls.paired_files to be a non-empty list of file paths.")
        else:
            # ensure each file exists
            for p in pf:
                if not p:
                    self.errors.append("Empty path found in paired control files list.")
                else:
                    if is_placeholder(p):
                        self.warnings.append(f"Paired control path looks like a placeholder: {p!r}")
                    if not os.path.exists(p):
                        # fatal per your request
                        self.errors.append(f"Paired control file not found: {p}")

    # 5) If prebuilt or both: require control prebuilt config (version or local dir)
    if method in ('prebuilt', 'both'):
        pre = (cfg.get('prebuilt') or {})
        ctrl = (pre.get('control') or {})
        kg1 = (pre.get('kg1') or {})

        ctrl_local = (ctrl.get('control_hash_local_dir') or "").strip()
        ctrl_version = (ctrl.get('control_hash_version') or "").strip()

        if not ctrl_local and not ctrl_version:
            self.errors.append("controls.method includes 'prebuilt' but prebuilt.control has neither control_hash_local_dir nor control_hash_version (one is required).")
        else:
            if ctrl_local:
                if is_placeholder(ctrl_local):
                    self.warnings.append("prebuilt.control.control_hash_local_dir looks like a placeholder; ensure you set a real path.")
                elif not os.path.exists(ctrl_local):
                    # fatal because user requested local prebuilt in effect (we're in prebuilt mode and local provided)
                    self.errors.append(f"prebuilt.control.control_hash_local_dir not found on host: {ctrl_local}")

            if ctrl_version and not isinstance(ctrl_version, str):
                self.warnings.append("prebuilt.control.control_hash_version should be a string like 'v1.0'")

        # KG1 handling: required unless user set no_kg1_removal = True
        if not cfg.get('no_kg1_removal', False):
            kg1_local = (kg1.get('kg1_hash_local_dir') or "").strip()
            kg1_version = (kg1.get('kg1_hash_version') or "").strip()
            if not kg1_local and not kg1_version:
                self.errors.append("KG1 removal enabled but prebuilt.kg1 has neither kg1_hash_local_dir nor kg1_hash_version (one is required) — or set no_kg1_removal: true to skip KG1 removal.")
            else:
                if kg1_local:
                    if is_placeholder(kg1_local):
                        self.warnings.append("prebuilt.kg1.kg1_hash_local_dir looks like a placeholder; ensure you set a real path.")
                    elif not os.path.exists(kg1_local):
                        self.errors.append(f"prebuilt.kg1.kg1_hash_local_dir not found on host: {kg1_local}")

    # 6) Mode-specific checks (region)
    mode = cfg.get('mode')
    is_region = (mode == 'region') or ('window' in cfg and cfg.get('window'))
    if is_region:
        window = cfg.get('window', {}) or {}
        wsize = window.get('window_size')
        rfile = window.get('region_file')
        if not wsize and not rfile:
            self.errors.append("Region mode requires window.window_size or window.region_file to be set.")
        if rfile:
            if is_placeholder(rfile):
                self.warnings.append("window.region_file looks like a placeholder; replace with a real regions file")
            elif not os.path.exists(rfile):
                # fatal per your instruction
                self.errors.append(f"window.region_file not found on host: {rfile}")

    # 7) Tools/resource-specific checks (Slurm and container)
    tools = (cfg.get('tools') or {})
    container_type = tools.get('container_type')
    resource_mgmt = tools.get('resource_mgmt')

    if resource_mgmt == 'slurm':
        slurm = (cfg.get('resources') or {}).get('slurm') or {}
        for key in ('account', 'partition', 'time'):
            val = slurm.get(key)
            if not val:
                self.errors.append(f"resources.slurm.{key} is required when resource_mgmt == 'slurm'")
            else:
                if is_placeholder(val):
                    self.warnings.append(f"resources.slurm.{key} looks like a placeholder; replace with your cluster-specific value: {key}={val!r}")

        # optional: enforce sing image existence if singularity will be used on the submit host (fatal)
        img = (cfg.get('container') or {}).get('singularity_image')
        if container_type == 'singularity':
            if not img:
                self.errors.append("container.singularity_image is required when tools.container_type == 'singularity'")
            else:
                if is_placeholder(img):
                    self.warnings.append("container.singularity_image appears to be a placeholder")
                elif not os.path.exists(img):
                    # fatal — if you require the image present on the host
                    self.errors.append(f"Singularity image not found on host: {img}")

    # 8) Container/docker checks when container_type == 'docker'
    if container_type == 'docker':
        img = (cfg.get('container') or {}).get('docker_image')
        if not img:
            self.errors.append("container.docker_image is required when tools.container_type == 'docker'")
        else:
            if is_placeholder(img):
                self.warnings.append("container.docker_image appears to be a placeholder")

    # 9) Other helpful checks
    # If method == 'prebuilt' but user also supplied paired files at top-level, we already required they exist.
    # If inputs.paired_control_files present but method == 'paired_files' it's already enforced above.
    # warn if no KG1 removal configured but user set no_kg1_removal true (just echo)
    if cfg.get('no_kg1_removal', False):
        self.warnings.append("no_kg1_removal is true: KG1 cohort removal will be skipped.")

    # end of validate: self.errors / self.warnings populated


    def plan(self):
        # run the full sequence
        self.select_profile()
        if self.errors:
            return
        self.merge()
        self.detect_overrides()
        self.resolve_prebuilt()
        # inject selected_profile into effective for reporting
        if self.selected_profile_name:
            self.effective['_selected_profile'] = self.selected_profile_name
        self.validate()

    # ---------- reporting utilities ----------
    def summary_text(self):
        out = []
        prof = self.selected_profile_name or "(none)"
        out.append(f"Selected profile: {prof}")
        method = (self.effective.get('controls') or {}).get('method')
        out.append(f"Resolved controls.method: {method}")
        # prebuilt summary
        pre = self.effective.get('_resolved', {}).get('prebuilt', {})
        for k in ('control', 'kg1'):
            r = pre.get(k, {})
            if not r:
                out.append(f"prebuilt.{k}: not configured")
                continue
            if r.get('use_local'):
                out.append(f"prebuilt.{k}: using LOCAL dir -> {r.get('local_dir')}")
            else:
                out.append(f"prebuilt.{k}: will FETCH version {r.get('fetch_version')}, recommended local path: {r.get('recommended_local_path')}")
        return "\n".join(out)

    def print_report(self, verbose=False):
        print_header("RUFUS PLAN REPORT")
        print(self.summary_text())
        if self.warnings:
            print_header("WARNINGS")
            for w in self.warnings:
                print(f"- {w}")
        if self.errors:
            print_header("ERRORS")
            for e in self.errors:
                print(f"- {e}")
        if verbose:
            print_header("EFFECTIVE CONFIG (YAML)")
            print(yaml.safe_dump(self.effective, sort_keys=False))

# ---------- CLI ----------
def parse_args():
    p = argparse.ArgumentParser(description="RUFUS planner: parse config, resolve profile, validate, and report.")
    p.add_argument("--config", "-c", required=True, help="Path to run-config.yaml")
    p.add_argument("--profile", "-p", required=False, help="Profile to use (overrides default_profile)")
    p.add_argument("--workdir", "-w", required=False, default=os.getcwd(), help="Working directory (planner uses to suggest fetch paths)")
    p.add_argument("--describe-profile", action="store_true", help="Print schema hints for the named profile and exit")
    p.add_argument("--generate-template", metavar="PROFILE", help="Emit a minimal template YAML for a profile")
    p.add_argument("--validate", action="store_true", help="Validate configuration; exit non-zero on errors")
    p.add_argument("--dump-effective", metavar="OUTFILE", help="Write the fully resolved (merged) config to a YAML file")
    p.add_argument("--verbose", action="store_true", help="Print verbose effective config in report")
    return p.parse_args()

def describe_profile(cfg, prof_name):
    profiles = (cfg.get('profiles') or {})
    if prof_name not in profiles:
        print(f"Profile '{prof_name}' not found. Available: {list(profiles.keys())}")
        return 2
    prof = profiles[prof_name] or {}
    print_header(f"PROFILE: {prof_name}")
    print("Fields present in this profile (key: example / note):")
    def walk(d, prefix=[]):
        for k, v in (d or {}).items():
            path = '.'.join(prefix + [k])
            print(f"- {path}: {v!r}")
            if isinstance(v, dict):
                walk(v, prefix + [k])
    walk(prof)
    print("\nHelpful hints (schema):")
    for key,h in SCHEMA_HINTS.items():
        if key.startswith(tuple(prof.keys())) or any(key.startswith(p+".") for p in prof.keys()):
            print(f"- {key}: {h[1]}")
    return 0

def generate_template(cfg, profile_name):
    # minimal template: include top-level required keys + the profile-specific keys
    base = {}
    # include a few top-level helpful defaults (from cfg if present)
    for key in ('inputs','algorithm','working_dir','threads_per_job','prebuilt','controls'):
        if key in cfg:
            base[key] = cfg[key]
    profiles = (cfg.get('profiles') or {})
    if profile_name not in profiles:
        print(f"# ERROR: profile '{profile_name}' not found", file=sys.stderr)
        sys.exit(2)
    prof = profiles[profile_name] or {}
    # overlay profile keys only (so template is focused)
    out = deep_merge(base, prof)
    # ensure controls.method appears even if null -> helpful to user
    if 'controls' not in out:
        out['controls'] = {'method': None}
    # print template YAML to stdout
    yaml.safe_dump(out, sys.stdout, sort_keys=False)

def main():
    args = parse_args()
    cfg = load_yaml(args.config)
    if args.describe_profile:
        if not args.profile:
            print("Please supply --profile <name> with --describe-profile", file=sys.stderr)
            sys.exit(2)
        sys.exit(describe_profile(cfg, args.profile))

    if args.generate_template:
        generate_template(cfg, args.generate_template)
        return

    planner = Planner(cfg, cli_profile=args.profile, workdir=args.workdir)
    planner.plan()
    planner.print_report(verbose=args.verbose)

    # write effective config if requested
    if args.dump_effective:
        with open(args.dump_effective, 'w') as fh:
            yaml.safe_dump(planner.effective, fh, sort_keys=False)
        print(f"Wrote resolved config to {args.dump_effective}")

    if planner.errors:
        print("\nPlanner found fatal errors; aborting (exit 3).", file=sys.stderr)
        sys.exit(3)
    if args.validate:
        if planner.warnings:
            print("\nValidation completed with warnings (exit 0).")
        else:
            print("\nValidation successful (no warnings).")
    # otherwise exit 0 for success
    return

if __name__ == "__main__":
    main()
