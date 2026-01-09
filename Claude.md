# KimeraCS - Claude.md

A Windows Forms application for viewing and editing 3D models from Final Fantasy VII.

## Tech Stack
- .NET 9.0 Windows Forms
- OpenTK 5.0.0-pre.15 (OpenGL 3.3+ with GLSL shaders)
- ~45,600 lines of C# code

### OpenGL Namespaces
The codebase uses two OpenTK namespaces:
- **`OpenTK.Graphics.OpenGL`** - Modern OpenGL 3.3+ (used by GLRenderer, Mesh, ShaderProgram)
- **`OpenTK.Graphics.OpenGL.Compatibility`** - Legacy/immediate mode GL (used by ModelDrawing, Utils, Lighting, FF7*.cs)

**Important:** Legacy GL functions in Compatibility namespace require suffixes:
- `GL.Translate()` → `GL.Translated()`
- `GL.Scale()` → `GL.Scaled()`
- `GL.Rotate()` → `GL.Rotated()`
- `GL.Color3()` → `GL.Color3f()`
- `GL.Color4()` → `GL.Color4f()`
- `GL.Vertex3()` → `GL.Vertex3f()`
- `GL.Normal3()` → `GL.Normal3f()`
- `GL.TexCoord2()` → `GL.TexCoord2f()`
- `GL.MultMatrix()` → `GL.MultMatrixd()`

## Build Commands
```bash
# Build from repository root
dotnet build src/KimeraCS.csproj

# Run the application
dotnet run --project src/KimeraCS.csproj
```

## Project Structure
```
src/
├── Core/                    # Core utilities and modern rendering
│   ├── GLRenderer.cs        # Modern OpenGL renderer with mesh caching
│   ├── Mesh.cs              # GPU mesh classes (PModelMesh, GroupMesh, LineMesh, etc.)
│   ├── ShaderProgram.cs     # Shader compilation and uniforms
│   ├── ModelDrawing.cs      # Vertex/polygon operations, GLRenderer wrappers
│   ├── MatrixStack.cs       # CPU-side matrix stack (replaces GL matrix stack)
│   └── VisualizationHelpers.cs  # Modern mesh builders for debug visuals
├── Shaders/                 # GLSL shaders (model, line, point .vert/.frag)
├── Properties/              # Settings and resources
├── resources/               # UI assets (icons, cursors)
├── frmSkeletonEditor.cs     # Main window - skeleton/animation viewing
├── frmPEditor.cs            # Polygon editor - vertex/texture/color editing
├── frm*.cs                  # Various dialog forms
└── FF7*.cs                  # File format handlers
```

## Key Namespaces
| Namespace | Purpose |
|-----------|---------|
| `KimeraCS` | Main UI forms and FF7 format handlers |
| `KimeraCS.Core` | Utilities and enums |
| `KimeraCS.Rendering` | Modern OpenGL rendering (GLRenderer, Mesh, ShaderProgram) - located in Core/ folder |

## FF7 File Formats

### Skeletons
| Format | Extension | Description |
|--------|-----------|-------------|
| HRC | `.HRC` | Field character skeletons (text-based, hierarchical) |
| AA/DA | `.D` | Battle character/enemy skeletons (binary) |
| RSD | `.RSD` | Resource references linking models to bones |

### Models
| Format | Extension | Handler | Description |
|--------|-----------|---------|-------------|
| P-Model | `.P` | `FF7PModel.cs` | Main FF7 3D model format |
| TMD | (varies) | `FF7TMDModel.cs` | PlayStation TMD format |
| 3DS | `.3DS` | `Model_3DS.cs` | Import from Autodesk 3D Studio |

### Animations
| Format | Extension | Handler | Description |
|--------|-----------|---------|-------------|
| A-File | `.A` | `FF7FieldAnimation.cs` | Field character animations |
| DA-File | `DA` suffix | `FF7BattleAnimation.cs` | Battle animation packs |
| A00 | `.A00` | Similar to DA | Magic model animations |

### Textures
| Format | Extension | Handler | Description |
|--------|-----------|---------|-------------|
| TEX | `.TEX` | `FF7TEXTexture.cs` | FF7 texture format (paletted or direct) |

## Core Classes

### Data Structures
| Class | File | Purpose |
|-------|------|---------|
| `FF7PModel` | FF7PModel.cs | P-format 3D models (vertices, polygons, groups, textures) |
| `FF7TMDModel` | FF7TMDModel.cs | TMD format models |
| `FF7Skeleton` | FF7Skeleton.cs | Static orchestrator for skeleton loading/rendering |
| `FF7FieldSkeleton` | FF7FieldSkeleton.cs | Field character skeletons (.HRC) |
| `FF7BattleSkeleton` | FF7BattleSkeleton.cs | Battle/enemy skeletons |
| `FF7FieldAnimation` | FF7FieldAnimation.cs | Field animations (.A files) |
| `FF7BattleAnimation` | FF7BattleAnimation.cs | Battle animations |
| `FF7TEXTexture` | FF7TEXTexture.cs | TEX texture format |

### Rendering
| Class | File | Purpose |
|-------|------|---------|
| `GLRenderer` | Core/GLRenderer.cs | Modern OpenGL renderer with mesh caching (vertex colors + polygon colors) |
| `Mesh`/`GroupMesh`/`PModelMesh` | Core/Mesh.cs | GPU mesh abstraction (VAO/VBO/EBO) |
| `LineMesh`/`PointMesh` | Core/Mesh.cs | Line and point rendering for debug visuals |
| `ShaderProgram` | Core/ShaderProgram.cs | Shader compilation and uniforms |
| `MatrixStack`/`MatrixManager` | Core/MatrixStack.cs | CPU-side matrix stack (replaces GL.PushMatrix/PopMatrix) |
| `VisualizationHelpers` | Core/VisualizationHelpers.cs | Creates meshes for normals, bounding boxes, axes, wireframes |
| `SkeletonRenderer` | Core/SkeletonRenderer.cs | Modern skeleton bone rendering (computes transforms on CPU) |
| `ModelDrawing` | Core/ModelDrawing.cs | Vertex/polygon operations, matrix syncing wrappers for GLRenderer |
| `Lighting` | Lighting.cs | 4 configurable lights; `SetLightsModern()` for shaders, `SetLights()` for legacy |

### UI Forms
| Form | Purpose |
|------|---------|
| `FrmSkeletonEditor` | Main editor - skeleton/animation viewing |
| `FrmPEditor` | Polygon editor with paint tools |
| `FrmTextureViewer` | UV coordinate editing |
| `FrmFieldDB` | Field character database browser |
| `FrmBattleDB` | Battle model database |
| `FrmMagicDB` | Magic model database |
| `FrmGroupProperties` | Group blend/shading properties |

### Utilities
| Class | File | Purpose |
|-------|------|---------|
| `Utils` | Utils.cs | Math (Point2D, Point3D, Quaternion), vectors, bounding boxes |
| `FileTools` | FileTools.cs | Config file handling, database management |
| `UndoRedo` | UndoRedo.cs | Skeleton editor undo/redo |
| `UndoRedoPE` | UndoRedoPE.cs | Polygon editor undo/redo |

## Configuration

**File:** `Kimera.cfg` in application root

Key settings:
- LGP paths (CHAR.LGP, BATTLE.LGP, MAGIC.LGP source/destination)
- Working folders for models, animations, textures
- Window positions/sizes
- Undo buffer capacity
- Default interpolation frames

## P Editor Display Modes

The polygon editor (`frmPEditor`) has three display modes controlled by `drawMode`. All modes use the modern GLRenderer:

| Mode | Value | Description |
|------|-------|-------------|
| K_MESH | 0 | Wireframe mesh only (black lines) |
| K_PCOLORS | 1 | Polygon colors + wireframe overlay |
| K_VCOLORS | 2 | Vertex colors with textures |

### GLRenderer Drawing Functions
| Function | Purpose |
|----------|---------|
| `DrawPModelModern` | Renders with vertex colors and textures |
| `DrawPModelWireframe` | Renders wireframe with solid override color |
| `DrawPModelPolygonColors` | Renders with polygon colors (flat shading, no textures) |

The polygon color mesh uses a separate cache (`PColorMeshCacheByName`) since it builds vertices differently (duplicates vertices per-polygon with the polygon's color).

### Matrix Setup Functions (Utils.cs)

| Function | Purpose |
|----------|---------|
| `SetCameraPModel` | Sets projection + modelview matrices, syncs to GLRenderer |
| `SetCameraModelView` | Sets modelview with translate/rotate/scale (calls LoadIdentity) |
| `ConcatenateCameraModelView` | Concatenates additional transforms onto existing modelview |
| `ConcatenateCameraModelViewQuat` | Like above but uses quaternion, calls LoadIdentity first |

**Important:** Matrix setup in `DrawPModelEditor` must match `MouseMove` to avoid visual misalignment between display modes. Both should use `SetCameraPModel(..., 1,1,1)` followed by `ConcatenateCameraModelView` for resize.

## Rendering Pipeline

### Blend Modes
| Mode | Value | Effect |
|------|-------|--------|
| BLEND_AVG | 0 | 50% background + 50% polygon |
| BLEND_ADD | 1 | Additive blending |
| BLEND_SUB | 2 | Subtractive blending |
| BLEND_25P | 3 | 25% polygon blending |
| BLEND_NONE | 4 | No blending (opaque) |

### Mesh Structure
```
PModel (data) → PModelMesh (GPU)
    └─ GroupMesh[] (per-group VAO/VBO/EBO)
        ├─ Position, Normal, UV, Color buffers
        └─ TextureID, BlendMode, AlphaRef metadata

PModelMesh.FromPModel(model, usePolygonColors):
  - usePolygonColors=false: Vertex colors from Model.Vcolors[]
  - usePolygonColors=true:  Polygon colors from Model.Pcolors[] (same color for all 3 verts)
```

## Key Patterns

- **Vertex indexing:** Polygons store relative vertex indices in `Polys[i].Verts[0..2]`. To get absolute vertex index: `Model.Polys[polyIdx].Verts[vertIdx] + Model.Groups[groupIdx].offsetVert`
- **Static classes** for global state (FF7PModel, FF7Skeleton, Utils, FileTools)
- **Structs** for binary serialization compatibility with FF7 formats
- **MarshalAs** attributes for fixed-size arrays in binary formats
- **Circular buffers** for undo/redo (default capacity: 10)
- **Per-context resources** in GLRenderer for multi-window support

## Common Tasks

### Loading a model
```
File → FF7PModel.ReadPModel() → Load textures → GLRenderer creates mesh → Display
```

### Editing workflow
```
Select item → UndoRedo snapshot → Modify → Refresh view → (Undo restores snapshot)
```

### Rendering a frame
```
GLRenderer.DrawPModelModern() → Get/create cached mesh → For each group: set blend, bind texture, draw
```

## OpenTK 5.x API Notes

### Enum Renames
| Old (OpenTK 4.x) | New (OpenTK 5.x) |
|------------------|------------------|
| `MaterialFace` | `TriangleFace` |
| `CullFaceMode` | `TriangleFace` |
| `EnableCap.Texture2D` | `EnableCap.Texture2d` |
| `TextureTarget.Texture2D` | `TextureTarget.Texture2d` |
| `BufferUsageHint` | `BufferUsage` |

### Modern GL API Changes (OpenTK.Graphics.OpenGL)
| Old | New |
|-----|-----|
| `GL.GetProgram(handle, param, out val)` | `val = GL.GetProgrami(handle, param)` |
| `GL.GetShader(handle, param, out val)` | `val = GL.GetShaderi(handle, param)` |
| `GL.Uniform1(loc, int)` | `GL.Uniform1i(loc, int)` |
| `GL.Uniform1(loc, float)` | `GL.Uniform1f(loc, float)` |
| `GL.UniformMatrix4(loc, transpose, ref mat)` | `GL.UniformMatrix4f(loc, count, transpose, ref mat)` |
| `GL.GenTextures(1, array)` | `texId = GL.GenTexture()` |

### Texture ID Types
OpenTK 5.x expects `int` for texture IDs. If storing as `uint`, cast when calling GL functions:
```csharp
GL.BindTexture(TextureTarget.Texture2d, (int)texId);
GL.IsTexture((int)texId);
```

### PixelFormat Ambiguity
When using Compatibility namespace with `System.Drawing.Imaging`, fully qualify OpenGL's PixelFormat:
```csharp
OpenTK.Graphics.OpenGL.Compatibility.PixelFormat.Bgra
```

## Rendering

### Model Shader (model.frag)
The fragment shader uses standard Lambertian diffuse lighting to match legacy OpenGL fixed-function behavior:
- **Lighting enabled**: `ambient + diffuse` where diffuse = `max(dot(normal, lightDir), 0.0)`
- **Lighting disabled**: Raw vertex/polygon color (no brightness boost)
- **Override color mode**: For wireframe rendering, bypasses all color/texture calculations

Key uniforms:
| Uniform | Type | Purpose |
|---------|------|---------|
| `useTexture` | bool | Sample texture or use white |
| `enableLighting` | bool | Apply lighting calculations |
| `useOverrideColor` | bool | Use solid color (wireframe mode) |
| `overrideColor` | vec3 | Color when override is enabled |
| `baseAlpha` | float | Alpha multiplier for blend modes |

### Multi-Light System (Modern)
The shader supports 4 lights (indices 0-3):
- **0 (Right)**: Positioned at +X, 50% intensity
- **1 (Left)**: Positioned at -X, 50% intensity
- **2 (Front)**: Positioned at +Z, 100% intensity (default)
- **3 (Rear)**: Positioned at -Z, 75% intensity

```csharp
GLRenderer.LightPositions[index] = new Vector3(x, y, z);
GLRenderer.LightColors[index] = new Vector3(r, g, b);
GLRenderer.LightEnabled[index] = true/false;
```

### Matrix Syncing (Hybrid Rendering)
When calling modern rendering from within legacy matrix contexts, you MUST sync the legacy GL matrices to GLRenderer:

```csharp
// Get legacy matrices
double[] projMatrix = new double[16];
GL.GetDouble(GetPName.ProjectionMatrix, projMatrix);
double[] mvMatrix = new double[16];
GL.GetDouble(GetPName.ModelviewMatrix, mvMatrix);

// Save and set GLRenderer matrices
var savedProjection = GLRenderer.ProjectionMatrix;
var savedView = GLRenderer.ViewMatrix;
var savedModel = GLRenderer.ModelMatrix;

GLRenderer.ProjectionMatrix = ToMatrix4(projMatrix);
GLRenderer.ViewMatrix = Matrix4.Identity;
GLRenderer.ModelMatrix = ToMatrix4(mvMatrix);

// Call modern rendering
GLRenderer.DrawLinesModern(mesh);

// Restore matrices
GLRenderer.ProjectionMatrix = savedProjection;
GLRenderer.ViewMatrix = savedView;
GLRenderer.ModelMatrix = savedModel;
```

**Functions that require this pattern:** `DrawBox`, `ShowNormals`, `SkeletonRenderer.RenderFieldSkeletonBones`, `SkeletonRenderer.RenderBattleSkeletonBones`
