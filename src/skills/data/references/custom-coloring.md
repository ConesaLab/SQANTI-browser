# Custom coloring

Override track colors with `--my-palette palette.json`. Missing file → exit 1.

```bash
python -m sqanti_browser --gtf corrected.gtf --classification classification.txt \
    --output my_hub --genome hg38 --my-palette my_colors.json
```

Colors apply everywhere: hub tracks, track description HTML pages, and the interactive table
reports (including SVG diagrams). Validation-track colors are embedded in the bigBed, so they
persist if you copy tracks into a UCSC custom track.

## JSON format

Three optional top-level keys. Each value is an RGB array `[R, G, B]` **or** a hex string
`"#RRGGBB"` (the `#` is optional).

- `standard` — base category colors
- `highlight` — color for the top-FL isoform per group
- `validation_tracks` — colors for `CAGE_peaks`, `polyA_peaks`, `star_sj`, `reference`

**Category keys** (use these exact strings under `standard`/`highlight`):
`full-splice_match`, `incomplete-splice_match`, `novel_in_catalog`, `novel_not_in_catalog`,
`genic`, `antisense`, `fusion`, `intergenic`, `genic_intron`.

## Partial specification

Everything is optional — omit any key and it keeps its default. If you provide `standard`
but omit `highlight`, highlight colors are **auto-derived by darkening** the standard colors.
To turn off highlighting entirely, pass `--no-highlight`.

## Minimal example

Only the listed categories change; the rest keep defaults, and highlight is auto-darkened.

```json
{
    "standard": {
        "full-splice_match": "#FF0000",
        "incomplete-splice_match": [0, 255, 0]
    }
}
```

## Full example (with validation tracks)

```json
{
    "standard": {
        "full-splice_match": [107, 174, 214],
        "novel_in_catalog": [120, 198, 121]
    },
    "highlight": {
        "full-splice_match": [69, 111, 137]
    },
    "validation_tracks": {
        "CAGE_peaks": [0, 128, 0],
        "polyA_peaks": [200, 0, 0],
        "star_sj": [21, 101, 192],
        "reference": [70, 70, 70]
    }
}
```

See `example/example_palette.json` for a complete reference file.

## Default colors

Standard (CSS name → hex): full-splice_match steelblue #6BAED6, incomplete-splice_match
coral #FC8D59, novel_in_catalog lightgreen #78C679, novel_not_in_catalog crimson #D62F4B,
genic darkgray #969696, antisense mediumaquamarine #66C2A4, fusion goldenrod #DAA520,
intergenic darksalmon #E9967A, genic_intron cadetblue #41B6C4.

Highlight defaults darken each: steelblue, peru, forestgreen, darkred, dimgray, teal,
darkgoldenrod, sienna, darkslategray.

Validation defaults (RGB): CAGE_peaks 200,133,41 (amber, TSS); polyA_peaks 17,76,145 (blue,
TTS); star_sj 207,228,205 (light green, splice junctions); reference 70,70,70 (gray).
