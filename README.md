# VL.Nurbsy (WIP)

A NURBS Bezier library for [vvvv](https://vvvv.org) and [Stride](https://github.com/stride3d/stride). Based on [LNLib](https://github.com/BIMCoderLiang/LNLib) by [BIMCoderLiang](https://github.com/BIMCoderLiang).

This project is a C# port and adaptation of [LNLib](https://github.com/BIMCoderLiang/LNLib) by [BIMCoderLiang](https://github.com/BIMCoderLiang), aims to bring [The NURBS Book 2nd Edition](https://link.springer.com/book/10.1007/978-3-642-97385-7) algorithms to the native C# domain.

⚠️ **This project is currently in active development** and may be subject to breaking changes, deprecations, and unstable features.

## Packages

- **Nurbsy** - Pure C# Core library containing all geometric algorithms. Can be used in any Stride Game or .NET application. 
- **VL.Nurbsy** - Wrapper for **vvvv gamma**. Provides a node set to visualize and manipulate curves and surfaces interactively. 

### Curves
- **NurbsCurve** - Evaluation and manipulation of non-uniform rational B-spline curves.
- **BezierCurve** - Specialized support for Bezier curves.

### Surfaces
- **NurbsSurface** - Evaluation and manipulation of non-uniform rational B-spline surfaces.
- **BezierSurface** - Specialized support for Bezier patches.

### Installation
```sh
// vvvv
nuget install VL.Nurbsy 

// Stride
dotnet add package Nurbsy
```

### Getting Started

Help patches for [vvvv](https://vvvv.org) are available via Help Browser.

## Contributing

Contributions are welcome! 
Open [issues](https://github.com/antokhio/VL.Nurbsy/issues) or submit pull requests.
Questions welcome on [vvvv forum](https://forum.vvvv.org).

## License
The source code is published under [LGPL 2.1](https://www.gnu.org/licenses/), the license is available [here](LICENSE).
