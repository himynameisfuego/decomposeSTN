# decomposeSTN
Fuzzy decomposition of sound into sines, transients and noise. Available both for Matlab and Python.

* L. Fierro, and V. Välimäki. _"**Enhanced Fuzzy Decomposition of Sound Into Sines, Transients and Noise.**"_.  Accepted for publication into the Journal of the Audio Engineering Society, 2023.

Link: [ArXiv](https://arxiv.org/abs/2210.14041)

## Abstract

The decomposition of sounds into sines, transients, and noise is a long-standing research problem in audio processing. The current solutions for this three-way separation detect either horizontal and vertical structures or anisotropy and orientations in the spectrogram to identify the properties of each spectral bin and classify it as sinusoidal, transient, or noise. This paper proposes an enhanced three-way decomposition method based on fuzzy logic, enabling soft masking while preserving the perfect reconstruction property. The proposed method allows each spectral bin to simultaneously belong to two classes, sine and noise or transient and noise. Results of a subjective listening test against three other techniques are reported, showing that the proposed decomposition yields a better or comparable quality. The main improvement appears in transient separation, which enjoys little or no loss of energy or leakage from the other components and performs well for test signals presenting strong transients.

## Dependencies

### Matlab
* [Signal Processing Toolbox](https://www.mathworks.com/products/signal.html)

### Python
* SciPy
* NumPy

## Usage
### Recommended parameters

``Fs = 44100 # Hz

nWin1 = 8192 # samples 

nWin2 = 512 # samples ``

### Matlab
[xs, xt, xn] = decomposeSTN(audioInput,Fs,[nWin1 nWin2]);

### Python
import decomposeSTN as STN

[xs, xt, xn] = STN.decSTN(audioInput,Fs,[nWin1 nWin2])

## Contributing
Suggestions and contributions to the code are both welcomed and encouraged. Please open an issue to discuss your changes and submit a pull request.

## License
Please refer to [**LICENCE.md**](LICENSE.md) for further information.
