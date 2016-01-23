library(ggplot2)

# MULTIPLE DATA SETS
#ori <- read.table('testing/rotate_oriloc_ecoli.txt', header=FALSE)
multi <- read.table('testing/rotate_multiscale_staph.txt', header=FALSE)
zcurve <- read.table('testing/rotate_orifinder_staph.txt', header=FALSE)
data <- rbind(cbind(multi, Program="Multiscaling"), cbind(zcurve, Program="Ori-Finder"))


pdf("rotate_multi-orifinder_staph.pdf", height=12, width=12)
ggplot(data) + geom_point(aes(x = V1, y = V2, color = Program), alpha = 0.5) + 
	xlab("Offset (bp)") + ylab("Distance From 0 (bp)") + 
	ggtitle("Overlay of Multiscaling, & Ori-Finder Rotation Tests with S. Aureus") + theme_bw() +
	scale_color_brewer(palette = "Set1")
dev.off()


# SINGLE DATA SET
#data <- read.table('testing/rotate_multiscale_staph.txt', header=FALSE)
#pdf("rotate_overlay_staph.png", height=12, width=12)
#ggplot(data) + geom_point(aes(x = V1, y = V2)) + xlab("Offset (bp)") + ylab("Distance From 0 (bp)") + ggtitle("Multiscaling Rotation Tests with S. Aureus") + theme_bw()
#dev.off()