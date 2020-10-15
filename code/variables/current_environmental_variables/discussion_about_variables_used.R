#Explanation about what variables are interesting to include


############################
##### PREDICTORES ##########
############################
###LA DECISION  DE SI O NO INCLUIAR EN LOS ANALSISI DE SELECCION DE VARIABLES ESTA EXPLICADA EN notebook_pine_distribution_models.tex.

##Variables por buscar. Intenta seleccionar variables que se sepa que tienen interés sobre la especie, no las elijas al tun tun. 
    #Clima:
        #SI: Obtuvimos el indice de humedad (precipiatacion-evaporacion potencial) del grupo de nick. Con ese indice y la tempreratura minima y maxima, creamos todas las variables de bioclim usando biovars de dismo package. Para más info mira notebook_pine_distribution_models.tex.
        #NO: Radiacion solar (kJ m-2 day-1): De http://www.worldclim.org/version2. WorldClim version 2 has average monthly climate data for minimum, mean, and maximum temperature and for precipitation for 1970-2000.
    #Topografía:
        #SI (SOLO RESAMPLING):elevacion (metros de altitud)
            #Vamos a usar los datos que tiene Nick en su intranet de GTOPO30, con una resolucion de 1 km en el ecuador. Hay que preguntarle a Nick si esa resolución es suficiente, o hace falta más. En Cuyo caso tendríamos que descargar los datos de SRTM y ASTER, con una resolucion de 30 metros. SRTM falta a partir de los 60 grados para arriba y abajo (respectivamente en cada hemisferio), pero eso se completa con ASTER. 
            #Los datos de STRTM están en el link http://dds.cr.usgs.gov/srtm/. Ojo que es la versión 2.1 con los huecos rellenados mediante interpolación. En esa página hay que seleccionar version2_1/SRTM30/ para ir a la versión 2 y datos con una resolucion de 30 metros. En ese directorio tienes una carpeta para cada zona muestreados. Dentro de la caperta de cada zona, te interesa el archivo DEM, dentor hay un archivo .DEM que puedes abrir con R usando la funcion raster. Una vez descargados y cargados todos, habrá que juntar todos los raster con la funcon merge de raster (http://stackoverflow.com/questions/15876591/merging-multiple-rasters-in-r).
        #NO: Topographic Position Index (diversidad topográfica?): TPI is the difference between the value of a cell and the mean value of its 8 surrounding cells.
        #NO: slope (pendiente)
    #NO: Huella humana (0-100) (http://sedac.ciesin.columbia.edu)
        #De una página de la NASA (http://sedac.ciesin.columbia.edu/data/set/wildareas-v2-human-footprint-geographic) te descargas un zip de los datos mundiales y en la carpeta hf_v2geo coges el archivo w001001x.adf, lo abres con la funcion raster y listo. OJO, tienes que dejar todos los archvios que vienen con él, sino no lo puedes leer. 
        #The Human Influence Index (HII) is a global dataset of 1-kilometer grid cells, created from nine global data layers covering human population pressure (population density), human land use and infrastructure (built-up areas, nighttime lights, land use/land cover), and human access (coastlines, roads, railroads, navigable rivers). The dataset is produced by the Wildlife Conservation Society (WCS) and the Columbia University Center for International Earth Science Information Network (CIESIN) and is available in the Geographic Coordinate system.
    #NO: Cover land: 
        ##La capa de cubierta vegetal puede ser descargada desde http://due.esrin.esa.int/page_globcover.php. Se llama globcover y te presenta diferentes CATEGORIAS de cubiertas. En la carpeta de los datos (Globcover2009_V2.3_Global_) hay un excel llamado Globcover2009_Legend.xls donde está la legenda de los colores. Se ve como cada tipo de cubierta tiene uan composicion RGB (red, green, blue) diferente. Si es un Rainfed croplands tiene mas rojo y verde que azul, si es Closed to open (>15%) broadleaved forest regularly flooded tiene mas verde y menos rojo. No tengo muy claro si usar esta varaible, habrá que hablar con Nick. 
    #Suelo: 
        #SI: SoilGrid. Soilgrid represents a collection of updatable soil property and class maps of the world at 1 km and 250 m spatial resolution produced using automated soil mapping based on machine learning algorithms. Segun Rapha, tiene puntos reales y luego rellena los huecos con interpolación espacial, igual que Bioclim. http://www.isric.org/content/soilgrids
            #Para saber el nombre de cada variable en el ftp tienes que entrar en https://soilgrids.org/#/?layer=geonode:taxnwrb_250m, buscar la capa de la varible que te interesa, y luego darle a download, ahí te saldrá la ID de la variable. 
            #Características del sitio: 
                #Absolute depth to bedrock (in cm). I think that depth rooting diferences between pines species could be important to cope with drought, but I don`t know if depth of soil in general can be important for pine species distribution. En ftp.soilgrids.org vas a data/aggregated/1km/ busca BDTICM_M_sl lo que sea.
            #Fisicas: Textura del suelo, es un factor que determina la distribucion, sobre todo en algunas especies como Pinus pinea. 
                #Clay content (%): En ftp.soilgrids.org vas a data/aggregated/1km/ para conseguir los datos de todo el planeta agregados, con una resolucion de 1 km. Clay content es CLYPPT_M_slX, donde X es un número del 1 al 7, siendo el dato a 0, 5, 15, 30, 60, 100 y 200 cm. 
                #Silt content (%). Contenido en limo. En ftp.soilgrids.org vas a data/aggregated/1km/ busca SLTPPT_M_sl lo que sea.
                #Sand content (%). Contenido en arena. En ftp.soilgrids.org vas a data/aggregated/1km/ busca SNDPPT_M_sl lo que sea. 
            #Quimicas.  
                #pH (index * 10): En ftp.soilgrids.org vas a data/aggregated/1km/ busca PHIHOX_M_sl lo que sea. Hay algunos pinos con ciertos requerimientos de pH.  
                #Capacidad de intercambio cationico (cmolc/kg). En En ftp.soilgrids.org vas a data/aggregated/1km/ busca CECSOL_M_sl... Puede afectar a la biodisponibilidad de ciertos nutrientes. https://es.wikipedia.org/wiki/Capacidad_de_intercambio_cati%C3%B3nico
                #Carbono orgánico (g/kg). En ftp.soilgrids.org vas a data/aggregated/1km/ busca ORCDRC_M_sl... No termino de verlo claro. 
