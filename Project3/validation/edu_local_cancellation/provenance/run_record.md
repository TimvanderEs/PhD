# Exact rerun record

- EC2 run root: `/home/ec2-user/HELIOS_validation/edu_local_cancellation_20260918/run_validated_20260918`
- Master process: PID 83976
- Started: 18 September 2026 at approximately 17:08 UTC
- Corrected run completed: 19 September 2026 at 13:25:35 UTC
- Historical-control run completed: 19 September 2026 at 14:02:52 UTC
- Post-processing root: `/home/ec2-user/HELIOS_validation/edu_local_cancellation_20260918/analysis_validated_20260919`
- Runtime: R 4.2.2 on Amazon Linux 2023
- LAVA: 0.1.5
- LAVA execution: one sequential locus worker per configuration; the control
  and corrected configurations ran concurrently.

The copied `.log` and `.sessionInfo.txt` files are stored beside the exact
control and corrected outputs under `02_lava_inputs/rerun/results`. Input and
reference hashes are in `exact_run_input_checksums.sha256`; copied output
hashes are in `exact_run_outputs.sha256`; post-processing hashes are in
`analysis_outputs.sha256`.
