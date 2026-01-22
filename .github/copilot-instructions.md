# Copilot Instructions

## General Guidelines
- First general instruction
- Second general instruction

## Code Style
- Use specific formatting rules
- Follow naming conventions

## Project-Specific Rules
- In the VL.Nurbsy project, use immutable collections with reference equality for change detection.
- The `NurbsSurfaceNode` base class handles `Build()` via the `Invalidate()/Update()` pattern called each frame.
- Conversion from `T` to `ControlPoint<T>` is intentional to bridge the user-facing simple API with internal weighted control points.