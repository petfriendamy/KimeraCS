# KimeraCS - Claude.md

A Windows Forms application for viewing and editing 3D models from Final Fantasy VII.

## Tech Stack
- .NET 9.0 Windows Forms
- OpenTK 4.8.2 (OpenGL 3.3+ with GLSL shaders)
- ~45,600 lines of C# code

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
├── Core/                    # Core utilities, ModelDrawing (vertex/polygon ops), BlendModes
├── Rendering/               # Modern OpenGL pipeline (GLRenderer, Mesh, ShaderProgram)
├── Shaders/                 # GLSL shaders (model.vert/frag, line.vert/frag)
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
| `KimeraCS.Rendering` | Modern OpenGL rendering (GLRenderer, Mesh classes) |

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
| `GLRenderer` | Rendering/GLRenderer.cs | Modern OpenGL renderer with mesh caching |
| `Mesh`/`GroupMesh`/`PModelMesh` | Rendering/Mesh.cs | GPU mesh abstraction (VAO/VBO/EBO) |
| `ShaderProgram` | Rendering/ShaderProgram.cs | Shader compilation and uniforms |
| `ModelDrawing` | Core/ModelDrawing.cs | Vertex/polygon operations (paint, select, pick), legacy GL compatibility |
| `Lighting` | Lighting.cs | 4 configurable lights |

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

The polygon editor (`frmPEditor`) has three display modes controlled by `drawMode`:

| Mode | Value | Renderer | Description |
|------|-------|----------|-------------|
| K_MESH | 0 | Legacy GL | Wireframe mesh only |
| K_PCOLORS | 1 | Legacy GL | Polygon colors + wireframe overlay |
| K_VCOLORS | 2 | GLRenderer | Vertex colors with textures (modern shaders) |

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
