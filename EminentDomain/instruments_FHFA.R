#### A LOT OF INPUT !!!!!!

z <- data.frame(FHFA$numpanels1x_noreligion, 
                FHFA$numpanels1x_jd_public, 
                FHFA$numpanels1x_dem, 
                FHFA$numpanels2x_dem, 
                FHFA$numpanels3x_dem, 
                (FHFA$numpanels1x_dem)^2, 
                (FHFA$numpanels1x_dem)^3, 
                FHFA$numpanels1x_female, 
                FHFA$numpanels2x_female,
                FHFA$numpanels3x_female, 
                FHFA$numpanels1x_nonwhite,
                FHFA$numpanels2x_nonwhite,
                FHFA$numpanels3x_nonwhite, 
                FHFA$numpanels1x_black,
                FHFA$numpanels2x_black,
                FHFA$numpanels3x_black, 
                FHFA$numpanels1x_jewish,
                FHFA$numpanels2x_jewish,
                FHFA$numpanels3x_jewish, 
                FHFA$numpanels1x_catholic,
                FHFA$numpanels2x_catholic,
                FHFA$numpanels3x_catholic, 
                FHFA$numpanels2x_noreligion,
                FHFA$numpanels3x_noreligion, 
                FHFA$numpanels1x_instate_ba,
                FHFA$numpanels2x_instate_ba,
                FHFA$numpanels3x_instate_ba, 
                FHFA$numpanels1x_ba_public,
                FHFA$numpanels2x_ba_public,
                FHFA$numpanels3x_ba_public, 
                FHFA$numpanels2x_jd_public,
                FHFA$numpanels3x_jd_public, 
                (FHFA$numpanels1x_jd_public)^2/5,
                (FHFA$numpanels1x_jd_public)^3/25, 
                FHFA$numpanels1x_elev,
                FHFA$numpanels2x_elev,
                FHFA$numpanels3x_elev, 
                (FHFA$numpanels1x_elev)^2,
                (FHFA$numpanels1x_elev)^3, 
                FHFA$numpanels1x_female_black,
                FHFA$numpanels2x_female_black,
                FHFA$numpanels3x_female_black, 
                FHFA$numpanels1x_female_noreligion,
                FHFA$numpanels2x_female_noreligion,
                FHFA$numpanels3x_female_noreligion, 
                FHFA$numpanels1x_black_noreligion,
                FHFA$numpanels2x_black_noreligion,
                FHFA$numpanels3x_black_noreligion, 
                FHFA$numpanels1x_female_jd_public,
                FHFA$numpanels2x_female_jd_public, 
                FHFA$numpanels3x_female_jd_public, 
                FHFA$numpanels1x_black_jd_public,
                FHFA$numpanels2x_black_jd_public, 
                FHFA$numpanels3x_black_jd_public, 
                FHFA$numpanels1x_noreligion_jd_pub,
                FHFA$numpanels2x_noreligion_jd_pub, 
                FHFA$numpanels3x_noreligion_jd_pub, 
                FHFA$numpanels1x_mainline,
                FHFA$numpanels2x_mainline,
                FHFA$numpanels3x_mainline, 
                FHFA$numpanels1x_evangelical,
                FHFA$numpanels2x_evangelical, 
                FHFA$numpanels3x_evangelical, 
                FHFA$numpanels1x_protestant,
                FHFA$numpanels2x_protestant,
                FHFA$numpanels3x_protestant, 
                FHFA$numpanels1x_dem * FHFA$numpanels1x_female,
                FHFA$numpanels1x_dem *  FHFA$numpanels1x_black, 
                FHFA$numpanels1x_dem * FHFA$numpanels1x_jewish,
                FHFA$numpanels1x_dem *  FHFA$numpanels1x_catholic, 
                FHFA$numpanels1x_dem * FHFA$numpanels1x_noreligion,
                FHFA$numpanels1x_dem *  FHFA$numpanels1x_instate_ba, 
                FHFA$numpanels1x_dem * FHFA$numpanels1x_jd_public,
                FHFA$numpanels1x_dem *  FHFA$numpanels1x_mainline, 
                FHFA$numpanels1x_dem * FHFA$numpanels1x_evangelical,
                FHFA$numpanels1x_dem *  FHFA$numpanels1x_protestant, 
                FHFA$numpanels1x_dem * FHFA$numpanels1x_ba_public,
                FHFA$numpanels1x_dem *  FHFA$numpanels1x_elev, 
                FHFA$numpanels1x_dem * FHFA$numpanels1x_nonwhite, 
                FHFA$numpanels1x_female *  FHFA$numpanels1x_jewish,
                FHFA$numpanels1x_female *  FHFA$numpanels1x_catholic, 
                FHFA$numpanels1x_female *  FHFA$numpanels1x_noreligion,
                FHFA$numpanels1x_female *  FHFA$numpanels1x_instate_ba, 
                FHFA$numpanels1x_female *  FHFA$numpanels1x_jd_public,
                FHFA$numpanels1x_female *  FHFA$numpanels1x_mainline, 
                FHFA$numpanels1x_female *  FHFA$numpanels1x_evangelical,
                FHFA$numpanels1x_female *  FHFA$numpanels1x_protestant, 
                FHFA$numpanels1x_female *  FHFA$numpanels1x_ba_public,
                FHFA$numpanels1x_female *  FHFA$numpanels1x_elev, 
                FHFA$numpanels1x_female *  FHFA$numpanels1x_nonwhite, 
                FHFA$numpanels1x_jd_public *  FHFA$numpanels1x_jewish,
                FHFA$numpanels1x_jd_public *  FHFA$numpanels1x_catholic, 
                FHFA$numpanels1x_jd_public *  FHFA$numpanels1x_noreligion,
                FHFA$numpanels1x_jd_public *  FHFA$numpanels1x_instate_ba, 
                FHFA$numpanels1x_jd_public *  FHFA$numpanels1x_mainline, 
                FHFA$numpanels1x_jd_public *  FHFA$numpanels1x_evangelical,
                FHFA$numpanels1x_jd_public *  FHFA$numpanels1x_protestant, 
                FHFA$numpanels1x_jd_public *  FHFA$numpanels1x_ba_public,
                FHFA$numpanels1x_jd_public *  FHFA$numpanels1x_elev, 
                FHFA$numpanels1x_jd_public *  FHFA$numpanels1x_nonwhite, 
                FHFA$numpanels1x_ba_public *  FHFA$numpanels1x_jewish,
                FHFA$numpanels1x_ba_public *  FHFA$numpanels1x_catholic, 
                FHFA$numpanels1x_ba_public *  FHFA$numpanels1x_noreligion,
                FHFA$numpanels1x_ba_public *  FHFA$numpanels1x_instate_ba, 
                FHFA$numpanels1x_ba_public *  FHFA$numpanels1x_mainline, 
                FHFA$numpanels1x_ba_public *  FHFA$numpanels1x_evangelical,
                FHFA$numpanels1x_ba_public *  FHFA$numpanels1x_protestant, 
                FHFA$numpanels1x_ba_public *  FHFA$numpanels1x_elev, 
                FHFA$numpanels1x_ba_public *  FHFA$numpanels1x_nonwhite, 
                FHFA$numpanels2x_jd_public *  FHFA$numpanels1x_jewish,
                FHFA$numpanels2x_jd_public *  FHFA$numpanels1x_catholic, 
                FHFA$numpanels2x_jd_public *  FHFA$numpanels1x_noreligion,
                FHFA$numpanels2x_jd_public *  FHFA$numpanels1x_instate_ba, 
                FHFA$numpanels2x_jd_public *  FHFA$numpanels1x_mainline, 
                FHFA$numpanels2x_jd_public *  FHFA$numpanels1x_evangelical,
                FHFA$numpanels2x_jd_public *  FHFA$numpanels1x_protestant, 
                FHFA$numpanels2x_jd_public *  FHFA$numpanels1x_ba_public,
                FHFA$numpanels2x_jd_public *  FHFA$numpanels1x_elev, 
                FHFA$numpanels2x_jd_public *  FHFA$numpanels1x_nonwhite, 
                FHFA$numpanels3x_jd_public *  FHFA$numpanels1x_jewish,
                FHFA$numpanels3x_jd_public *  FHFA$numpanels1x_catholic, 
                FHFA$numpanels3x_jd_public *  FHFA$numpanels1x_noreligion, 
                FHFA$numpanels3x_jd_public *  FHFA$numpanels1x_instate_ba, 
                FHFA$numpanels3x_jd_public *  FHFA$numpanels1x_mainline, 
                FHFA$numpanels3x_jd_public *  FHFA$numpanels1x_evangelical,
                FHFA$numpanels3x_jd_public *  FHFA$numpanels1x_protestant, 
                FHFA$numpanels3x_jd_public *  FHFA$numpanels1x_ba_public,
                FHFA$numpanels3x_jd_public *  FHFA$numpanels1x_elev, 
                FHFA$numpanels3x_jd_public *  FHFA$numpanels1x_nonwhite, 
                FHFA$numpanels2x_female *  FHFA$numpanels1x_jewish,
                FHFA$numpanels2x_female *  FHFA$numpanels1x_catholic, 
                FHFA$numpanels2x_female *  FHFA$numpanels1x_noreligion,
                FHFA$numpanels2x_female *  FHFA$numpanels1x_instate_ba, 
                FHFA$numpanels2x_female *  FHFA$numpanels1x_jd_public,
                FHFA$numpanels2x_female *  FHFA$numpanels1x_mainline, 
                FHFA$numpanels2x_female *  FHFA$numpanels1x_evangelical,
                FHFA$numpanels2x_female *  FHFA$numpanels1x_protestant, 
                FHFA$numpanels2x_female *  FHFA$numpanels1x_ba_public,
                FHFA$numpanels2x_female *  FHFA$numpanels1x_elev, 
                FHFA$numpanels2x_female *  FHFA$numpanels1x_nonwhite, 
                FHFA$numpanels2x_ba_public *  FHFA$numpanels1x_jewish,
                FHFA$numpanels2x_ba_public *  FHFA$numpanels1x_catholic, 
                FHFA$numpanels2x_ba_public *  FHFA$numpanels1x_noreligion,
                FHFA$numpanels2x_ba_public *  FHFA$numpanels1x_instate_ba, 
                FHFA$numpanels2x_ba_public *  FHFA$numpanels1x_mainline, 
                FHFA$numpanels2x_ba_public *  FHFA$numpanels1x_evangelical,
                FHFA$numpanels2x_ba_public *  FHFA$numpanels1x_protestant, 
                FHFA$numpanels2x_ba_public *  FHFA$numpanels1x_elev, 
                FHFA$numpanels2x_ba_public *  FHFA$numpanels1x_nonwhite, 
                FHFA$numpanels1x_instate_ba *  FHFA$numpanels1x_jewish,
                FHFA$numpanels1x_instate_ba *  FHFA$numpanels1x_catholic, 
                FHFA$numpanels1x_instate_ba *  FHFA$numpanels1x_noreligion,
                FHFA$numpanels1x_instate_ba *  FHFA$numpanels1x_instate_ba, 
                FHFA$numpanels1x_instate_ba *  FHFA$numpanels1x_mainline, 
                FHFA$numpanels1x_instate_ba *  FHFA$numpanels1x_evangelical,
                FHFA$numpanels1x_instate_ba *  FHFA$numpanels1x_protestant, 
                FHFA$numpanels1x_instate_ba *  FHFA$numpanels1x_elev, 
                FHFA$numpanels1x_instate_ba *  FHFA$numpanels1x_nonwhite, 
                FHFA$numpanels1x_elev *  FHFA$numpanels1x_jewish,
                FHFA$numpanels1x_elev *  FHFA$numpanels1x_catholic, 
                FHFA$numpanels1x_elev *  FHFA$numpanels1x_noreligion, 
                FHFA$numpanels1x_elev *  FHFA$numpanels1x_mainline, 
                FHFA$numpanels1x_elev *  FHFA$numpanels1x_evangelical,
                FHFA$numpanels1x_elev *  FHFA$numpanels1x_protestant, 
                FHFA$numpanels1x_elev *  FHFA$numpanels1x_nonwhite, 
                FHFA$numpanels1x_protestant *  FHFA$numpanels1x_jewish,
                FHFA$numpanels1x_protestant *  FHFA$numpanels1x_catholic, 
                FHFA$numpanels1x_protestant *  FHFA$numpanels1x_noreligion, 
                FHFA$numpanels1x_protestant *  FHFA$numpanels1x_mainline, 
                FHFA$numpanels1x_protestant *  FHFA$numpanels1x_evangelical,
                FHFA$numpanels1x_protestant *  FHFA$numpanels1x_protestant, 
                FHFA$numpanels1x_protestant *  FHFA$numpanels1x_nonwhite, 
                FHFA$numpanels1x_mainline *  FHFA$numpanels1x_jewish,
                FHFA$numpanels1x_mainline *  FHFA$numpanels1x_catholic, 
                FHFA$numpanels1x_mainline *  FHFA$numpanels1x_noreligion,
                FHFA$numpanels1x_mainline *  FHFA$numpanels1x_evangelical, 
                FHFA$numpanels1x_mainline *  FHFA$numpanels1x_mainline, 
                FHFA$numpanels1x_mainline *  FHFA$numpanels1x_nonwhite )

names(z) <- cbind('1x_noreligion',
                  '1x_jd_public',
                  '1x_dem',
                  '2x_dem',
                  '3x_dem', 
                  '1x_dem^2',
                  '1x_dem^3', 
                  '1x_female',
                  '2x_female',
                  '3x_female', 
                  '1x_nonwhite',
                  '2x_nonwhite',
                  '3x_nonwhite', 
                  '1x_black',
                  '2x_black',
                  '3x_black', 
                  '1x_jewish',
                  '2x_jewish',
                  '3x_jewish', 
                  '1x_catholic',
                  '2x_catholic',
                  '3x_catholic', 
                  '2x_noreligion',
                  '3x_noreligion', 
                  '1x_instate_ba',
                  '2x_instate_ba',
                  '3x_instate_ba', 
                  '1x_ba_public',
                  '2x_ba_public',
                  '3x_ba_public', 
                  '2x_jd_public',
                  '3x_jd_public', 
                  '1x_jd_public^2',
                  '1x_jd_public^3', 
                  '1x_elev',
                  '2x_elev',
                  '3x_elev', 
                  '1x_elev^2',
                  '1x_elev^3', 
                  '1x_female_black',
                  '2x_female_black',
                  '3x_female_black', 
                  '1x_female_noreligion',
                  '2x_female_noreligion',
                  '3x_female_noreligion', 
                  '1x_black_noreligion',
                  '2x_black_noreligion',
                  '3x_black_noreligion', 
                  '1x_female_jd_public',
                  '2x_female_jd_public',
                  '3x_female_jd_public', 
                  '1x_black_jd_public',
                  '2x_black_jd_public',
                  '3x_black_jd_public', 
                  '1x_noreligion_jd_pub',
                  '2x_noreligion_jd_pub',
                  '3x_noreligion_jd_pub', 
                  '1x_mainline',
                  '2x_mainline',
                  '3x_mainline', 
                  '1x_evangelical',
                  '2x_evangelical',
                  '3x_evangelical', 
                  '1x_protestant',
                  '2x_protestant',
                  '3x_protestant', 
                  '1x_dem.*1x_female',
                  '1x_dem.*1x_black', 
                  '1x_dem.*1x_jewish',
                  '1x_dem.*1x_catholic', 
                  '1x_dem.*1x_noreligion',
                  '1x_dem.*1x_ba_public', 
                  '1x_dem.*1x_jd_public',
                  '1x_dem.*1x_mainline', 
                  '1x_dem.*1x_evangelical',
                  '1x_dem.*1x_protestant', 
                  '1x_dem.*1x_ba_public',
                  '1x_dem.*1x_elev',
                  '1x_dem.*1x_nonwhite',
                  '1x_female.*1x_jewish',
                  '1x_female.*1x_catholic', 
                  '1x_female.*1x_noreligion',
                  '1x_female.*1x_ba_public', 
                  '1x_female.*1x_jd_public',
                  '1x_female.*1x_mainline', 
                  '1x_female.*1x_evangelical',
                  '1x_female.*1x_protestant', 
                  '1x_female.*1x_ba_public',
                  '1x_female.*1x_elev',
                  '1x_female.*1x_nonwhite',
                  '1x_jd_public.*1x_jewish',
                  '1x_jd_public.*1x_catholic', 
                  '1x_jd_public.*1x_noreligion',
                  '1x_jd_public.*1x_instate_ba', 
                  '1x_jd_public.*1x_mainline', 
                  '1x_jd_public.*1x_evangelical',
                  '1x_jd_public.*1x_protestant', 
                  '1x_jd_public.*1x_ba_public',
                  '1x_jd_public.*1x_elev', 
                  '1x_jd_public.*1x_nonwhite', 
                  '1x_ba_public.*1x_jewish',
                  '1x_ba_public.*1x_catholic', 
                  '1x_ba_public.*1x_noreligion',
                  '1x_ba_public.*1x_instate_ba', 
                  '1x_ba_public.*1x_mainline', 
                  '1x_ba_public.*1x_evangelical',
                  '1x_ba_public.*1x_protestant', 
                  '1x_ba_public.*1x_elev', 
                  '1x_ba_public.*1x_nonwhite', 
                  '2x_jd_public.*1x_jewish',
                  '2x_jd_public.*1x_catholic', 
                  '2x_jd_public.*1x_noreligion',
                  '2x_jd_public.*1x_instate_ba', 
                  '2x_jd_public.*1x_mainline', 
                  '2x_jd_public.*1x_evangelical',
                  '2x_jd_public.*1x_protestant', 
                  '2x_jd_public.*1x_ba_public',
                  '2x_jd_public.*1x_elev', 
                  '2x_jd_public.*1x_nonwhite', 
                  '3x_jd_public.*1x_jewish',
                  '3x_jd_public.*1x_catholic', 
                  '3x_jd_public.*1x_noreligion',
                  '3x_jd_public.*1x_instate_ba', 
                  '3x_jd_public.*1x_mainline', 
                  '3x_jd_public.*1x_evangelical',
                  '3x_jd_public.*1x_protestant', 
                  '3x_jd_public.*1x_ba_public',
                  '3x_jd_public.*1x_elev', 
                  '3x_jd_public.*1x_nonwhite', 
                  '2x_female.*1x_jewish',
                  '2x_female.*1x_catholic', 
                  '2x_female.*1x_noreligion',
                  '2x_female.*1x_ba_public', 
                  '2x_female.*1x_jd_public',
                  '2x_female.*1x_mainline', 
                  '2x_female.*1x_evangelical',
                  '2x_female.*1x_protestant', 
                  '2x_female.*1x_ba_public',
                  '2x_female.*1x_elev',
                  '2x_female.*1x_nonwhite',
                  '2x_ba_public.*1x_jewish',
                  '2x_ba_public.*1x_catholic', 
                  '2x_ba_public.*1x_noreligion',
                  '2x_ba_public.*1x_instate_ba', 
                  '2x_ba_public.*1x_mainline', 
                  '2x_ba_public.*1x_evangelical',
                  '2x_ba_public.*1x_protestant', 
                  '2x_ba_public.*1x_elev', 
                  '2x_ba_public.*1x_nonwhite', 
                  '1x_instate_ba.*1x_jewish',
                  '1x_instate_ba.*1x_catholic', 
                  '1x_instate_ba.*1x_noreligion',
                  '1x_instate_ba.*1x_instate_ba', 
                  '1x_instate_ba.*1x_mainline', 
                  '1x_instate_ba.*1x_evangelical',
                  '1x_instate_ba.*1x_protestant', 
                  '1x_instate_ba.*1x_elev', 
                  '1x_instate_ba.*1x_nonwhite',
                  '1x_elev.*1x_jewish',
                  '1x_elev.*1x_catholic', 
                  '1x_elev.*1x_noreligion',
                  '1x_elev.*1x_mainline', 
                  '1x_elev.*1x_evangelical',
                  '1x_elev.*1x_protestant', 
                  '1x_elev.*1x_nonwhite', 
                  '1x_protestant.*1x_jewish',
                  '1x_protestant.*1x_catholic', 
                  '1x_protestant.*1x_noreligion',
                  '1x_protestant.*1x_mainline', 
                  '1x_protestant.*1x_evangelical',
                  '1x_protestant.*1x_protestant', 
                  '1x_protestant.*1x_nonwhite', 
                  '1x_mainline.*1x_jewish',
                  '1x_mainline.*1x_catholic', 
                  '1x_mainline.*1x_noreligion', 
                  '1x_mainline.*1x_evangelical',
                  '1x_mainline.*1x_mainline', 
                  '1x_mainline.*1x_nonwhite')
