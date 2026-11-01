# Frozen production calibration

`mrscloud_press_te30_quantitative/` is a byte-identical copy of the seven Phase 6a
bundle files. Identity and scalar are unchanged. Ordinary Git is appropriate:
the original bundle is about 10 MiB; its largest file is about 6 MiB.

`provenance_sources.zip` adds immutable construction-source snapshots: 408 files,
36,427,564 uncompressed bytes, 13,885,482 compressed bytes. Each member retains its
original repository-relative path and is checked against the existing normalization
metadata's SHA256. No archive member is imported or executed. This preserves
source validation without requiring optional research scripts or ignored captures
in a production checkout. The original seven files and their metadata are never
rewritten. The added archive is under the same existing upstream source attribution
and provenance as the bundled MRSCloud/FID-A code.

The original reference_manifest.json contains one historical author-specific
absolute comparison path. It is inert provenance, never a runtime resource lookup.
It is deliberately retained because editing frozen metadata would change its hash.
This is a documented exception to the general no-personal-path commit policy.
Six archived upstream MRSCloud/FID-A text sources also contain historical example
paths/comments. Their bytes must remain intact for the original source hashes;
they are never executed. No proprietary LCModel material is included.

Do not regenerate, normalize, patch or independently rescale these responses.
See docs/quantitative_simulator_user_guide.md for default lookup and overrides.
