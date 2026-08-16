## Security Policy for eskaks

### Reporting a Vulnerability

If you discover a security vulnerability in `eskaks`, we would appreciate your help in resolving it. Please report vulnerabilities to us privately, through the repository's security tab, rather than by opening a public issue. We ask that you do not publicly disclose the details of any vulnerabilities until we have had the chance to address them.

A private report gives us the chance to ship a fix before the details are searchable. A public issue does the opposite: it hands out a working reproducer to everyone before anyone has a patch, which is exactly the situation this policy exists to avoid.

### Supported Versions

`eskaks` is at version `0.1.0` and no release has been published yet. That means there is currently nothing to backport to: the only supported code is the default branch, and security fixes land there directly. Please make sure you are building from an up to date checkout before reporting, because a bug you hit in an old working copy may already be fixed upstream.

Once tagged releases exist, this section will follow the same rule the rest of our tools use:

- **Latest version** (actively maintained)
- **One prior version** (security fixes when feasible)

We state the support window explicitly instead of leaving it implied, because "we support the latest version" means very different things to a user sitting on a six month old build, and they deserve to know before they file a report that we cannot act on.

### Security Updates

We are committed to ensuring the security of `eskaks`. Critical vulnerabilities will be fixed as soon as possible, and minor vulnerabilities will be fixed in a timely manner. For updates and security notifications, please follow the repository's release page.

### Security Best Practices

- **Keep dependencies up to date.** `eskaks` is pure Rust and shells out to nothing, so its entire attack surface is its own parsers plus the crates in `Cargo.lock`. `cargo update` and an advisory check against that lockfile are therefore the highest value security action a downstream packager can take.

- **Treat every input file as untrusted.** `eskaks` reads FASTA, VCF, GFF3 and gzipped variants of those, and in practice those files arrive from public archives, collaborators and pipeline steps nobody in this repository controls. Assume a malformed or hostile file can reach the parsers, which is precisely why the `fuzz/` targets exist. If you find an input that makes `eskaks` panic, loop forever or allocate without bound, that is worth reporting even if it looks like "just" a crash.

- **Do not run `eskaks` with elevated privileges unless strictly necessary.** The tool needs read access to the input files and write access to the output path, and nothing more. Running it as root buys you no functionality and turns any parser bug into a much worse day.

- **Think before publishing an HTML report.** The interactive report is self contained by design: it makes no network requests and pulls in no remote assets, so viewing one cannot leak your data to a third party. The flip side is that it embeds your results, including sequence and gene identifiers taken straight from the input, into the file itself. A report built from samples with identifiable names is as sensitive as those names are, so treat it accordingly before attaching it to an issue or a public page.

### Contact
- Visit our [repository security tab](https://github.com/PathoGenOmics-Lab/eskaks/tree/main?tab=security-ov-file) for more information on our security posture.
