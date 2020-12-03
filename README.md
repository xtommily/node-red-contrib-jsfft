# JS FFT
A Node-Red FFT and InvFFT implementaton.
Package made with code from https://github.com/dntj/jsfft.

## FFT

#### Input
msg.payload should contain an array of *n* signal samples
``` js
[0.52,0.41,0.27,0.06,-0.13,-0.34,-0.42,-0.51,-0.54,-0.49...]
```

#### Output
The outcome is an array of *n* objects, and depends on the chosen algorithm type. For FFT it is
``` js
[{"freq":0,"real":-0.43663841485977173,"imag":0},
{"freq":31.25,"real":-0.7358816862106323,"imag":-0.25845131278038025},
{"freq":62.5,"real":1.0576090812683105,"imag":0.7042925357818604}..]
```

and for Magnitude+Phase
``` js
[{"freq":0,"magnitude":0.22450639307498932,"phase":0},
{"freq":31.25,"magnitude":0.48251003926265185,"phase":24.83929511444238},
{"freq":62.5,"magnitude":1.440573933658928,"phase":-54.4667332040494}..]
```
## InvFFT
#### Input
msg.payload should contain an array of *n* corresponding objects
``` js
[{"freq":0,"real":-0.43663841485977173,"imag":0},
{"freq":31.25,"real":-0.7358816862106323,"imag":-0.25845131278038025},
{"freq":62.5,"real":1.0576090812683105,"imag":0.7042925357818604}..]
```
#### Output
Outgoing msg.payload will contain an array of *n* signal samples
``` js
[0.52,0.41,0.27,0.06,-0.13,-0.34,-0.42,-0.51,-0.54,-0.49..]
```
