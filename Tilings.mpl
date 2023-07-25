Tilings := module()

export GetTiling := proc(
    tilingtype::string,
    coloring::{list, string}:=FAIL,
    scale::[realcons,realcons] := [1.0,1.0],
    tilingstyle::identical("line", "polygon", "polygonoutline"):="polygonoutline",
    tilingoptions::list := [],
    extraplotoptions::list := [],
    $)
local obj, out;

    obj := TilingObject( tilingtype, tilingoptions );
    out := RenderTiling( obj, coloring, scale, tilingstyle, extraplotoptions);

    return out;

end proc; # GetTiling


export TilingObject := proc( tilingtype, tilingoptions, $ );

    if tilingtype = "square" then
        return SquareTiling(tilingoptions[]);
    elif tilingtype = "triangle" then
        return TriangleTiling(tilingoptions[]);
    elif tilingtype = "hexagon" then
        return HexTiling(tilingoptions[]);
    end if;

    error "tiling type %1 not supported", tilingtype;

end proc; # TilingObject

export RenderTiling := proc(
    obj,
    coloring,
    scale,
    tilingstyle,
    extraplotoptions, $)
local colors, polys, out, i, l;

    if type(coloring, 'string') then
        if ColorTools:-KnownPalette(coloring) then
            colors := [seq(
                ColorTools:-Color(cat(coloring, " ", i)),
                    i=1..obj:-colors) ];
        else
            colors := [ColorTools:-Color(coloring) $ obj:-colors];
        end if;
    elif type(coloring, 'list') then
        if nops(coloring) < obj:-colors or
            irem( nops(coloring), obj:-colors ) <> 0
        then
            error "need %2 colors for tiling", obj:-colors;
        end if;
        colors := map(ColorTools:-Color, coloring);
    else
        colors := map(ColorTools:-Color, plots:-setcolors()[1..obj:-colors]);
    end if;

    out := seq(
        plottools:-translate(
        plottools:-rotate(
        plottools:-scale(
            plottools:-polygon(obj:-shapes[l[1]], 'style'=tilingstyle,
                'color'=colors[l[2]])
        , scale[]), l[4] ), (scale *~ l[3])[])
            , l in obj:-locations);

    return plots:-display( out, 'style'=tilingstyle, 'view'=[0..8,0..8], 'scaling'='constrained', extraplotoptions[] );

end proc; # RenderTiling


local SquareTiling := proc(w::posint:=8,l::posint:=8)
local i,j;

    return Record(shapes=[[[0, 0], [1, 0], [1, 1], [0, 1]]], colors=2,
        locations=[ # shape, color, pos, rot
                    seq(seq( [1, (i+j mod 2)+1, [i,j], 0],
                        i=0..w-1), j=0..l-1) ] );

end proc;


local TriangleTiling := proc(w::posint:=11,l::posint:=6)
local i,j;

    return Record(shapes=[[[0, 3/2], [-1/2*3^(1/2), 0], [1/2*3^(1/2), 0]]],
            'colors'=2,
            locations=[ # shape, color, pos, rot
                        seq(seq( [1, (i+j mod 2)+1, [i*1/2*3^(1/2),(3/2)*j+(i+j mod 2)*(3/2)], (i+j mod 2)*Pi],
                            i=0..w-1), j=0..l-1) ] );

end proc;


local HexTiling := proc(w::posint:=6,l::posint:=6)
local i, j;

    return Record(
            shapes=[ [[1, 3^(1/2)], [0, 3^(1/2)], [-1/2, 1/2*3^(1/2)], [0, 0], [1, 0], [3/2, 1/2*3^( 1/2)]] ],
            'colors'=3,
            locations=[ # shape, color, pos, rot
                seq(
                seq( [1, ((i mod 2)+j mod 3)+1, [i*3/2, j*sqrt(3) -(i mod 2)*sqrt(3)/2], 0], j=0..l-1),
                i=0..l-1)
            ] );

end proc;

end module:
