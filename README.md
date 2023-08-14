# AP_scripts_peterslab
Scripts written after opening lab

## Pipelines: 

#### Loading recordings
`ap.load_recording` - calls: 
- `ap.load_timelite/bonsai/mousecam/widefield/ephys`, based on what's available and workspace structure `load_parts`
- parses wheel with `ap.parse_wheel` into `wheel velocity` and binary `wheel_move` vectors

#### Ephys preprocessing
`ap.preprocess_neuropixels(animal,day)` - saves metadata, runs PyKilosort, translates spike times into Open Ephys timestamps

#### Data exploration
`ap.expscroll(wf_U,wf_V,wf_times,mousecam_fn,mousecam_times)` - scroll through widefield and/or mousecam
