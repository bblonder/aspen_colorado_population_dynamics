Attached are change grids for 14 species, representing the change between the reference period (1961-1990) and the decade around 2060.  

HOW THE MODELS WERE BUILT: Bioclimate models are based on species presence-absence data, historic climate, and slope/aspect data for over 850,000 points, including points sampled in FSveg vegetation polygons in southwestern Colorado, vegetation data from BLM, Mesa Verde National Park, Southern Ute Indian Tribe, and FIA plots in southwestern Colorado.  FIA plots across the southwestern US were selected and added to represent future climates that historically were not present in the study area.  These training data were fed into an algorithm called random forests to build multiple "forests" of many decision "trees".

VOTES:  Given climate and topographic data, all the "trees" vote as to whether the site is suitable for the species.  When less than 50% of trees vote for suitability in a cell in a given climate, it is considered unsuitable, otherwise suitable.

PROJECTIONS:  Projections are based on 3 general circulation models and 3 RCPs (carbon scenarios), resulting in 9 climates for the future.  The bioclimate model for each species was run with the reference climate and the 9 future climates for each grid cell.  The model output (votes for suitability) were averaged among the 9 future climate scenarios.

GRIDS:  The grids are geotiffs in lat/long WGS84 (CRS = +init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0).

Cells with the following values meet the given conditions:
0 neither suitable in the reference period nor in future.
1 Lost - suitable in reference period but very unsuitable in future (votes <30%).
2 Threatened - suitable in reference period but marginally unsuitable in future (votes >=30% and <50%).
3 Persistent - suitable in both periods.
4 Emergent - not suitable in reference period but suitable in future.
5 Water
6 Rock/Barren

Water and Rock/Barren were used to exclude species habitat where they occur.  Water was derived from USFS vegetation shapefiles.  Rock/Barren was derived from USFS and MVNP shapefiles.

WHAT IS SUITABILITY?  Bear in mind that the projections use climatic suitability, given the topography.  They do not take into account soil suitability (such as the lack thereof in alpine zones).  Strictly speaking, they do not project presence and absence.  Although historic suitability is closely related to species presence and absence, they could be less related in the future, due to lingering presence of species as the climate deteriorates or the failure of the species to migrate to emergent zones.

MAPS:  A png is supplied with an example of a projected map of change zones made from a change grid.

ACCURACY and UNCERTAINTY:  There are limits to the accuracy of the methods employed, and to the certainty attached to projections of the future.  

While the bioclimate models are quite effective at replicating the distribution of species at large scales, at smaller scales errors of omission (predicting absence where the species occurs) and commission (predicting presence where the species is absent) can be seen.  For examples, errors of omission may occur when a species is present in a small, suitable microsite in a grid cell that has, on average, unsuitable climate.  Although the best techniques were used, interpolating or downscaling climate information cannot replicate actual geographic variation in climate with complete accuracy.  

Although the projections are suitable for both regional and project-level planning, one should not place much reliance when planning at very small scales.  Projections for individual pixels should be considered less reliable than those for a larger area.

Climate projections are based on representative carbon pathways (RCP, i.e., emissions scenario) that may not represent the actual future trend in greenhouse gases.  For example, conditions projected for 2060 may actually occur sooner if emissions are higher than projected by RCP-6.0, or later if they are lower.  For these reasons, although boundaries between change zones must be precise for mapping and planning purposes, ideally they should be regarded as the best estimate of fuzzy boundaries, and the timing of the projected changes as likely but uncertain.

Please call or email if you have any questions or problems.

Jim Worrall			Suzanne Marchetti
jworrall@fs.fed.us		sbmarchetti@fs.fed.us
970-642-4453			970-642-4446

