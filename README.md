# STN decomposition
Fuzzy decomposition of sound into sines, transients and noise.

A MATLAB app for sines-transients-noise decomposition of audio signals. Developed using App Designer in Matlab 2020b.

* L. Fierro, and V. Välimäki. _"**Enhanced Fuzzy Decomposition of Sound Into Sines, Transients and Noise.**"_.  Proceedings of the 24th International Conference on Digital Audio Effects (DAFx20in21), Vienna, Austria.

## Abstract

Decomposition of sounds into their sinusoidal, transient, and noise components is an active research topic and a widely-used tool in audio processing. Multiple solutions have been proposed in recent years, using time-frequency representations to identify either horizontal and vertical structures or orientations and anisotropy in the spectrogram of the sound. 
This is SiTraNo: an easy-to-use MATLAB application with a graphic user interface for audio decomposition that enables visualization and access to the sinusoidal, transient, and noise classes, individually. This application allows the user to choose between different well-known separation methods to analyze an input sound file, to instantaneously control and remix its spectral components, and to visually check the quality of the decomposition, before producing the desired output file. The visualization of common artifacts, such as birdies and dropouts, is easy to get in SiTraNo. 

This app wants to promote experimenting with the sound decomposition process by observing the effect of variations for each spectral component on the original sound and by comparing different methods against each other, evaluating the separation quality both audibly and visually.

## Dependencies
### Matlab
* [Signal Processing Toolbox](https://www.mathworks.com/products/signal.html)

## Installation and use
* If your version of MATLAB is 2020b or later, download the latest [release](https://github.com/himynameisfuego/SiTraNo/releases/latest). If your version is 2020a or previous, refer to this [hotfix](https://github.com/himynameisfuego/SiTraNo/files/6351972/SiTraNo_HotFix_1.0.0.1.zip) until the next release.
* In MATLAB, navigate to the SiTraNo folder, open **SiTraNo.mlappinstall** and install. You will find SiTraNo in the "Apps" tab, in the "My apps" group. Click on it to execute the app.
* Upon launching SiTraNo, a navigation folder should pop up, asking you to choose the input audio file.

## Featured decomposition methods

* **HP**: Harmonic-Percussive separation [1]. Modes: hard mask, soft mask.
* **HPR**: Harmonic-Percussive-Residual separation [2]. Modes: single decomposition, two-round decomposition. 
* **ST**: Structure-Tensor-based separation [3].
* **Fuzzy**: Fuzzy logic decomposition [4].

## Contributing
Suggestions and contributions to the code are both welcomed and encouraged. Please open an issue to discuss your changes and submit a pull request.

## License
Please refer to [**LICENCE.md**](LICENSE.md) for further information.
