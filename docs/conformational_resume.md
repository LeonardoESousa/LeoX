# Resuming conformational searches

Run the same command again, or select option 6 and choose the same output
directory. The Gaussian template and xTB/CREST solvent must match the previous
search. Changing the Gaussian job concurrency limit is allowed.

Resume assumes that no previously submitted SLURM job remains active.

- xTB: reuse a valid optimized geometry only when the log records optimization
  convergence and normal termination. Otherwise discard incomplete outputs and restart in the same xTB folder.
- CREST: reuse a valid clustered ensemble only when the log records normal
  termination. Otherwise discard incomplete outputs and restart in the same CREST folder.
  Since conformer numbering can change, Gaussian files from that incomplete
  search are deleted.
- Gaussian: preserve existing inputs and completed opt/freq logs. Submit only
  unfinished initial jobs, then check completed results for a frequency retry.
- Frequency retries: overwrite the original Geometry-N-.com with
  Opt=ReadFC Guess=Read Geom=AllCheck followed by a separate frequency step.
  Temporarily retain the first log and checkpoint under Geometries/Retry/Originals
  for recovery and fallback; delete these recovery files when retries finish.
  Resume an interrupted second attempt using the original frequency Hessian.
  A completed second attempt is not repeated just because convergence is still
  imperfect. Incomplete or imaginary-frequency results are excluded as before.
- Classification: reuse reports when their hashes, source log hashes and
  thresholds match recorded completion metadata. Missing or changed reports
  trigger classification again. Older searches without this metadata are
  classified once to establish it.

Interrupted outputs are deleted or overwritten. Older Inputs/Interrupted and
Geometries/Interrupted archive folders are deleted on resume. Watcher always
retains checkpoints, including outside the conformational workflow. Progress is
appended to the existing workflow log.

Only conformation.csv, conformers_manifest.csv and conformers_unique.xyz are
written directly in the Conformational output folder. Resume state remains in
the stage folders. --classify-only retains its existing report-only behavior.
