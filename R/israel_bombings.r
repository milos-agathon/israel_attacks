# REPLICATION FILE - MILOS POPOVIC 5/13/2021
# MAKE A TIMELAPSE MAP OF BOMBINGS IN ISRAEL (2001-2019)!
# load packages
library(rgeos, quietly = T)
library(tweenr, quietly = T)
library(ggplot2, quietly=T) 
library(reshape2, quietly=T) 
library(grid, quietly=T) 
library(rgeos, quietly=T)
library(plyr, quietly=T)
library(dplyr, quietly=T)
library(rgdal, quietly=T)
library(animation, quietly=T)
library(scales, quietly=T)
library(data.table, quietly=T)
library(httr, quietly=T)
library(rmapshaper, quietly=T)
library(tidyr, quietly=T)
library(raster, quietly=T)
library(classInt, quietly=T)

set.seed(20210513)

# download and unzip files
url_cpg <- "https://github.com/justinelliotmeyers/official_israel_census_boundary_shapefile/raw/master/israel_demog2012.cpg.zip"
download.file(url_cpg, basename(url_cpg), mode="wb")
unzip("israel_demog2012.cpg.zip")

url_dbf <- "https://github.com/justinelliotmeyers/official_israel_census_boundary_shapefile/raw/master/israel_demog2012.dbf.zip"
download.file(url_dbf, basename(url_dbf), mode="wb")
unzip("israel_demog2012.dbf.zip")

url_prj <- "https://github.com/justinelliotmeyers/official_israel_census_boundary_shapefile/raw/master/israel_demog2012.prj.zip"
download.file(url_prj, basename(url_prj), mode="wb")
unzip("israel_demog2012.prj.zip")

url_sbn <- "https://github.com/justinelliotmeyers/official_israel_census_boundary_shapefile/raw/master/israel_demog2012.sbn.zip"
download.file(url_sbn, basename(url_sbn), mode="wb")
unzip("israel_demog2012.sbn.zip")

url_sbx <- "https://github.com/justinelliotmeyers/official_israel_census_boundary_shapefile/raw/master/israel_demog2012.sbx.zip"
download.file(url_sbx, basename(url_sbx), mode="wb")
unzip("israel_demog2012.sbx.zip")

url_shx <- "https://github.com/justinelliotmeyers/official_israel_census_boundary_shapefile/raw/master/israel_demog2012.shx.zip"
download.file(url_shx, basename(url_shx), mode="wb")
unzip("israel_demog2012.shx.zip")

url_shp <- "https://github.com/justinelliotmeyers/official_israel_census_boundary_shapefile/raw/master/israel_demog2012.shp.zip"
download.file(url_shp, basename(url_shp), mode="wb")
unzip("israel_demog2012.shp.zip")

# load shapefiles
isr <- readOGR(getwd(),
  "israel_demog2012", 
  verbose = TRUE, 
  stringsAsFactors = FALSE)

# calculate area size to be used for attack density
isr$area_sqkm <- area(isr) / 1000000

# load the data
# unfortunately, START GTD has no API service so 
# you have to download # the file, save it csv 
# format and import in R 
# below is the link to the most recent GTD file
df <- read.csv(file="https://www.dropbox.com/s/h037wmd8zf61gtq/globalterrorismdb_0221dist.csv?dl=1", 
			   header = T) %>%
	dplyr::filter(country_txt=="Israel" & iyear>=2001) %>% # only Israel
	dplyr::select(iyear, latitude, longitude) %>%
	na.omit() 

# let's determine to which settlement polygon every point belongs
isr$id <- 1:max(nrow(isr)) #first create a unique id for every polygon unit
d <- df[,c("longitude", "latitude")] # extract coordinate from df
coordinates(d) <- ~longitude+latitude
proj4string(d) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
i <- sp::over(d, isr) %>% # where the points fall
	cbind(df$iyear) # add iyear from df
names(i)[24] <- "iyear" # rename to iyear
i1 <- ddply(i, c("id", "iyear"), summarise, freq=length(id)) %>%
	na.omit()   # group by settlement and calculate sum of attacks
is <- merge(i1, isr, by="id", all=T) %>% # merge
      as.data.frame()

# fill in missing years for every location
yr_range <- seq(min(is$iyear, na.rm=T), 
         		max(is$iyear, na.rm=T), by=1)
yr_expanded <- expand.grid(yr_range, unique(isr$id))
colnames(yr_expanded) <- c("iyear", "id")
f <- merge(is, yr_expanded, by=c("id", "iyear"), all.y = T)

ff <- f %>%
group_by(id) %>%
complete(iyear = seq(min(iyear), 
         		max(iyear), by=1)) %>%
fill(id) %>%
fill(freq) %>%
ungroup() %>%
dplyr::select(id, iyear, freq, area_sqkm)

# we need to split the data.frame into list of years for tweenr
# in order to split the data, every unit must have every year
# otherwise the function returns an error
x <- split(ff, ff$iyear)

# use the list of years to set up transition features, number of frames etc.
tw <- tween_states(x, tweenlength= 2, statelength=3, 
                   ease='cubic-in-out',
                   nframes=60)
# attack density: attacks per square km of land area
tw$attack_dens <- tw$freq / tw$area_sqkm



# build a discrete legend
ni = classIntervals(tw$attack_dens, 
				   n = 5, 
				   style = 'quantile')$brks
# this function uses above intervals to create categories
labels <- c()
for(i in 1:length(ni)){
    labels <- c(labels, paste0(round(ni[i], 0), 
                             "–", 
                             round(ni[i + 1], 0)))
}
labels <- labels[1:length(labels)-1]

# finally, carve out the categorical variable 
# based on the breaks and labels above
tw$cat <- cut(tw$attack_dens, 
              breaks = ni, 
              labels = labels, 
              include.lowest = T)
levels(tw$cat) # let's check how many levels it has (5)

# label NAs, too
lvl <- levels(tw$cat)
lvl[length(lvl) + 1] <- "No attacks"
tw$cat <- factor(tw$cat, levels = lvl)
tw$cat[is.na(tw$cat)] <- "No attacks"
tw$cat <- revalue(tw$cat, c("0–0"="No attacks"))
levels(tw$cat)

#prepare the shapefile in data.frame format
israel <- fortify(isr, region = "id") %>% 
  mutate(id= as.numeric(id))

# i like to freeze the last frame for some time
# this code assigns 0.15 secs to each frame while
# extending the last to 4 secs
times <- c(rep(0.15, max(tw$.frame)-1), 4)
tw$iyear <- round(tw$iyear, 0) # remember to round years

# it's time to make our time-lapse map
oopt = ani.options(interval = .15)
saveGIF({for (i in 1:max(tw$.frame)) { # we need to define a loop to iterate through every frame
isr1 <- isr[,c(23)] # we don't want to duplicate columns to id is sufficient
isr1@data <- join(isr1@data, tw, by = "id") %>% # merge into new shapefile...
  filter(.frame==i) # ...and filter by frame
  map <- # map file to be created for every frame
    ggplot() +
  geom_map(data = israel, map = israel, # this is the base layer without colors
             aes(
              map_id = id),
             color = NA, 
             size = 0, 
             fill = NA) +
    geom_map(data = isr1@data, # this is the layer to be filled with values
             map = israel,
             aes(fill = cat, 
             map_id = id),
             color = NA, 
             size=0)  +
  scale_fill_manual(name= expression(paste("Attacks per ", km^{2}, "of land area")),
  values=c("grey80", '#ffd670', '#fb9164', '#dc5060', '#9a2d5c', '#422454'),
  labels=c("No attacks", "0-1", "1-2", "2-5", ">5"),
  drop=F)+
guides(fill=guide_legend(
            direction = "horizontal",
            keyheight = unit(2.5, units = "mm"),
            keywidth = unit(15, units = "mm"),
            title.position = 'top',
            title.hjust = 0.5,
            label.hjust = .5,
            nrow =1,
            byrow = T,
            #labels = labs,
            # also the guide needs to be reversed
            reverse = F,
            label.position = "bottom"
          )
    ) +
coord_map() +
expand_limits(x=israel$long,y=israel$lat)+
  labs(y="", 
        subtitle=paste0(as.character(as.factor(tail(tw %>% filter(.frame==i),1)$iyear))),
        title="Terrorist attacks in Israel (2001-2019)",
        caption="©2021 Milos Popovic (https://milospopovic.net)\n Data: GTD Database https://start.umd.edu/gtd/")+
theme_minimal() +
theme(plot.background = element_rect(fill = "white", color = NA), 
panel.background = element_rect(fill = "white", color = NA), 
legend.background = element_rect(fill = "white", color = NA),
legend.position = c(1, .05),
panel.border = element_blank(),
panel.grid.minor = element_blank(),
panel.grid.major = element_line(color = "white", size = 0),
plot.title = element_text(size=24, color="#422454", hjust=0.1, vjust=-15),
plot.subtitle = element_text(size=50, color="#f00065", hjust=0, vjust=-10, face="bold"),
plot.caption = element_text(size=14, color="grey60", hjust=.1, vjust=0),
legend.text = element_text(size=16, color="grey20"),
legend.title = element_text(size=22, color="grey20"),
strip.text = element_text(size=12),
plot.margin=unit(c(t=0, r=0, b=0, l=0), "cm"), # use this to cut extra blank space
axis.title.x = element_blank(),
axis.title.y = element_blank(),
axis.ticks = element_blank(),
axis.text.x = element_blank(),
axis.text.y = element_blank())
print(map)
  print(paste(i,"out of",max(tw$.frame))) # print progress
  ani.pause()}
},movie.name="israel_attacks.gif", 
dpi=600, 
ani.height=920, 
ani.width=920,
interval = times,
other.opts = "-framerate 10  -i image%03d.png -s:v 920x920 -c:v libx264 -profile:v high -crf 20  -pix_fmt yuv420p")