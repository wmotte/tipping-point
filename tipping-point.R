#Load libraries

library( 'ggplot2' )
library( 'reshape2' )
library( 'stringr' )
library( 'Matrix' )
library( 'plyr' )

library( 'gratia' )	# derivative functions
library( 'ggplot2' )	# plotting
library( "mgcv")        # GAM

#####################
# FUNCTIONS
#####################

# plotting function
number_ticks <- function( n ) { function( limits ) pretty( limits, n ) } 
row.selection <- NULL

#Read in static functioal connectivity data
read.static.data.z <- function()
{
	all.static.data <- NULL
	l <- NULL
	for(i in 1:10) #subjects
	{
	indir2 <- paste("directory with static FC folders for each subject", i, sep="")	
	matname <- paste0( indir2, '/static_FC_matrix.csv')
	print(matname)
	matrix <- read.csv(matname, row.names=1)
	l$rat <- i
	l$matrix <- matrix

	all.static.data <- rbind(all.static.data, l)
	}
	
	all.static.data <- as.data.frame(all.static.data, row.names=FALSE)
	return(all.static.data)
}

#### SCRIPTS #####	

#Set indir and outdir
indir <- ""

outdir <- ""
dir.create( outdir, showWarnings = FALSE )

#Read in functional and structural data
functional.static <- read.static.data.z()

for( i in 1:10) # For each individual dataset
{

average.functional <- data.frame(functional.static[i,2] )

structural <- read.csv(paste0(indir, "/structural_matrix.csv"), row.names=1)
rois <- read.table(paste0(indir, "/le_ri_HarvardOxford_no_id.txt")


#changes colnames and rownames
rownames(average.functional) <- colnames(average.functional) <-rownames(structural) <- colnames(structural) <- rois$V1

#Make lower triangle NA --> matrix is symmetrical (to avoid taking all connections twice)
average.functional[lower.tri(average.functional)] <- NA
structural[lower.tri(structural)] <- NA

#Melt datasets
melt.functional <- melt(as.matrix(average.functional), na.rm=TRUE)
melt.structural <- melt(as.matrix(structural), na.rm=TRUE)

#change colnames
colnames(melt.functional) <- c("ROI.1", "ROI.2", "FC")
colnames(melt.structural) <- c("ROI.1", "ROI.2", "SC")

#Calculate log-transform structural connectivity
melt.structural$log.SC <- log(melt.structural$SC)
melt.structural$log.SC[melt.structural$log.SC=="-Inf"] <- 0

#Combine datasets
data.all <- merge(melt.functional, melt.structural, by=c("ROI.1", "ROI.2"))

#Only select connections with an existing structural connectivity strength!
data.all <- data.all[data.all$SC > 0, ]

# GAM model fit --> fits non-linearity in the data
m <- gam( FC ~ s( log.SC, k = 5 ), data = data.all, method = "REML" )

## parameters for testing
n <- 500               # number of newdata values
EPS <- 1e-07           # finite difference
N <- 10000             # number of posterior draws (for CI estimation)

## where are we going to predict at? --> SC-values for which you are going to predict FC values
newd <- with( data.all, data.frame( log.SC = seq( min( log.SC ), max( log.SC ), length = n ) ) )

# fitted line --> predict values of FC with these SC values, based on the GAM model fit.
df_fit <- data.frame( x = newd$log.SC, y = predict( m, newdata = newd, type = 'response' ) )

# first derivatives of fitted GAM function
fd <- fderiv( m, newdata = newd, eps = EPS, unconditional = FALSE )

# get simultaneous confidence interval
set.seed( 42 )		# set the seed to make this repeatable 
sint <- confint( fd, type = "simultaneous", nsim = N )

# merge CI data with x
df_merged <- cbind( sint, x = newd$log.SC )

#Determine where the 95% CI does not contain 0 --> df_merged$zero == 0.; if 95% CI does contain 0 --> df_merged$zero == 1. 
interval <- cbind(df_merged$lower, df_merged$upper)
trueval <- 0
df_merged$zero <- apply(interval, 1, findInterval, x=trueval)

#Select rows where df_merged$zero changes from 1 to 0 or from 0 to 1 (potential tipping points)
rownumbers <- which(c(FALSE, tail(df_merged$zero, -1) != head(df_merged$zero ,-1)))

#Select the tipping point used for analyses: 
# - Do not use tipping points were df_merged == 1 --> This is when a 95% CI starts to contain 0 again 
# - When there is more than 1 tipping point; take the last tipping point
# - When there is no tipping point; just take the first row (x=0) --> delete those points later. 

selection <- NULL
if(length(rownumbers) > 1) {
for(r in 1:length(rownumbers))
	{
	row <- rownumbers[r]
	select <- df_merged[row,]
	select$rat <- i
	select$variable <- "static FC" 
	select$window <- "no" 
	selection <- rbind(selection, select)
	}
	
	selection2 <- selection[!(selection$zero==1),] 
	number <- nrow(selection2)
	select.final <- selection2[number,]
	
	row.selection <- rbind(row.selection, select.final)
	
} else if ( length(rownumbers)==1){
	row <- rownumbers[1]
	select <- df_merged[row,]
	if(select$zero == 0){
	select$rat <- i
	select$variable <- "static FC" 
	select$window <- "no" 
	} else if (select$zero ==1) {
	select <- df_merged[1,]
	select$rat <- i
	select$variable <- "static FC" 
	select$window <- "no" 
	} else {
	print("we have a problem!")
	}
	row.selection <- rbind(row.selection, select)
} else {
	select <- df_merged[1,]
	select$rat <- i
	select$variable <- "static FC" 
	select$window <- "no" 
	row.selection <- rbind(row.selection, select)
}

#Determine the y-value of the tipping point
select.yvalue <- df_merged[df_merged$x == select$x,]

# get plot: With tippoing point when there is a tipping point
if(select.yvalue$x > 0){
plot_1st <- ggplot( df_merged, aes( x = x, y = est ) ) +
	
	# original data
	geom_point( data = data.all, aes( x = log.SC, y = FC ), colour = "royalblue1", size = 1.5, alpha = 0.4 ) +
	# fitted line
	geom_line( data = df_fit, aes( x = x, y = y ), colour = "royalblue4", size = 2, alpha = 0.8 ) +
	# 1st derivative 
    geom_ribbon( aes( ymin = lower, ymax = upper ), alpha = 0.95, fill = '#fee6ce' ) +
	geom_line( size = 1.2, alpha = 0.4, colour = '#e6550d' ) +
	geom_hline( yintercept = 0, linetype = "dashed", color = "gray30", size = 0.5 ) +
	#Add tipping point as vertical line
	geom_segment(x=select$x, y=-0.7, xend=select$x, yend=select.yvalue$est, size = 0.8, linetype = "dashed") +
	# design
    labs( y = "Functional connectivity", x = "Structural connectivity" ) +
    scale_x_continuous( breaks = number_ticks( 10 ) ) + 
    scale_y_continuous( breaks = number_ticks( 10 ) ) +
	scale_fill_manual( values = c( "#1f78b4", "#E69F00", "gray50" ) ) +
    scale_color_manual( values = c( "#1f78b4", "#E69F00", "gray50" ) ) +
	theme_classic( base_size = 32 ) + 
	theme( legend.position = 'top', axis.title = element_text( face = "bold" ) )

} else {
plot_1st <- ggplot( df_merged, aes( x = x, y = est ) ) +
	
	# original data
	geom_point( data = data.all, aes( x = log.SC, y = FC ), colour = "royalblue1", size = 1.5, alpha = 0.4 ) +
	# fitted line
	geom_line( data = df_fit, aes( x = x, y = y ), colour = "royalblue4", size = 2, alpha = 0.8 ) +
	# 1st derivative 
    geom_ribbon( aes( ymin = lower, ymax = upper ), alpha = 0.95, fill = '#fee6ce' ) +
	geom_line( size = 1.2, alpha = 0.4, colour = '#e6550d' ) +
	geom_hline( yintercept = 0, linetype = "dashed", color = "gray30", size = 0.5 ) +
	# design
    labs( y = "Functional connectivity", x = "Structural connectivity" ) +
    scale_x_continuous( breaks = number_ticks( 10 ) ) + 
    scale_y_continuous( breaks = number_ticks( 10 ) ) +
	scale_fill_manual( values = c( "#1f78b4", "#E69F00", "gray50" ) ) +
    scale_color_manual( values = c( "#1f78b4", "#E69F00", "gray50" ) ) +
	theme_classic( base_size = 32 ) + 
	theme( legend.position = 'top', axis.title = element_text( face = "bold" ) )
}

plot_1st

# save plots to file
ggsave( file = paste0(  outdir, '/plot__1st_derivate_GAM_human', i, '.png' ), plot = plot_1st, dpi = 600, height = 8, width = 8 )
}


#Tipping points are the row.selection$x column
row.selection$tipping <- row.selection$x

#Remove rows where tipping point = 0 --> no tipping point
row.selection[row.selection$x==0, 'tipping' ] <- NA

#Save tipping points
write.csv(row.selection, file=paste0(outdir, '/tipping_points_human_whole_brain_static_FC.csv'))
