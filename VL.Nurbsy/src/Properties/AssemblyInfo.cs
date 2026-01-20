using Nurbsy;
using VL.Core.Import;

[assembly: ImportAsIs(Namespace = "VL", Category = "Nurbsy")]

[assembly: ImportType(typeof(CurveCurveIntersectionType), Category = "Nurbsy.Enums")]
[assembly: ImportType(typeof(LinePlaneIntersectionType), Category = "Nurbsy.Enums")]
[assembly: ImportType(typeof(CurveNormal), Category = "Nurbsy.Enums")]
[assembly: ImportType(typeof(SurfaceDirection), Category = "Nurbsy.Enums")]
[assembly: ImportType(typeof(SurfaceCurvature), Category = "Nurbsy.Enums")]
[assembly: ImportType(typeof(IntegratorType), Category = "Nurbsy.Enums")]
[assembly: ImportType(typeof(OffsetType), Category = "Nurbsy.Enums")]
