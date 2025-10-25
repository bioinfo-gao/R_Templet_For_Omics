# https://r-graph-gallery.com/210-custom-barplot-layout.html
# create dummy data
data <- data.frame(
  name=letters[1:5],
  value=sample(seq(4,15),5)
)

# Increase margin size
par(mar=c(11,4,4,4))

# Uniform color
barplot(height=data$value,
        col="#69b3a2",
        names.arg=c("very long group name 1",
                    "very long group name 2",
                    "very long group name 3",
                    "very long group name 4",
                    "very long group name 5"), 
        las=2 # means that the axis labels will be perpendicular to the axis
)

# las=0: Labels are always parallel to the axis (default).
# las=1: Labels are always horizontal.
# las=2: Labels are always perpendicular to the axis.
# las=3: Labels are always vertical.