:orphan:

****************
line properties
****************

  .. code-block:: C++

    color_map.insert(pair<string, string>("r", "'red'"));
    color_map.insert(pair<string, string>("g", "'green'"));
    color_map.insert(pair<string, string>("y", "'yellow'"));
    color_map.insert(pair<string, string>("b", "'blue'"));
    color_map.insert(pair<string, string>("w", "'white'"));
    color_map.insert(pair<string, string>("k", "'black'"));
    color_map.insert(pair<string, string>("c", "'cyan'"));


    style_map.insert(pair<string, string>("lp", "linespoints"));
    style_map.insert(pair<string, string>("l", "lines"));
    style_map.insert(pair<string, string>("p", "points"));
    style_map.insert(pair<string, string>("s", "steps"));
    style_map.insert(pair<string, string>("d", "dots"));
    style_map.insert(pair<string, string>("ip", "impulses"));
    style_map.insert(pair<string, string>("hs", "histeps"));
    style_map.insert(pair<string, string>("fs", "fsteps"));
    style_map.insert(pair<string, string>("fc", "filledcurves"));
    style_map.insert(pair<string, string>("hg", "histograms"));
    style_map.insert(pair<string, string>("box", "boxes"));

    line_marker_map.insert(pair<string, string>("+", "1"));
    line_marker_map.insert(pair<string, string>("x", "2"));
    line_marker_map.insert(pair<string, string>("*", "3"));
    line_marker_map.insert(pair<string, string>("@-", "4"));
    line_marker_map.insert(pair<string, string>("@", "5"));
    line_marker_map.insert(pair<string, string>("o-", "6"));
    line_marker_map.insert(pair<string, string>("o", "7"));
    line_marker_map.insert(pair<string, string>("^-", "8"));
    line_marker_map.insert(pair<string, string>("^", "9"));
    line_marker_map.insert(pair<string, string>("v-", "10"));
    line_marker_map.insert(pair<string, string>("v", "11"));
    line_marker_map.insert(pair<string, string>("#-", "12"));
    line_marker_map.insert(pair<string, string>("#", "13"));

    line_dash_map.insert(pair<string, string>(".", "'.'"));
    line_dash_map.insert(pair<string, string>("-", "'-'"));
    line_dash_map.insert(pair<string, string>("._", "'._'"));
    line_dash_map.insert(pair<string, string>("..-", "'..-'"));