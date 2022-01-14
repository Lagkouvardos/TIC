library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

input_dir = args[1]
output_pdf = args[2]

input_vector = paste(input_dir, '/alignment_vector.csv', sep='')

data_0 <- t(read.table(input_vector, header=FALSE, sep=','))

positions <- c(1:50000)
plot_data <- data.frame(positions, sums=data_0)

p <- ggplot(plot_data, aes(positions, sums))
g1 <- p + geom_point()
g3 <- g1 + xlab("SINA positions") + ylab("Counts")
g4 <- g3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
g5 <- g4 + labs(title = "Sum of SINA aligned positions")
g6 <- g5 + scale_x_continuous(breaks = scales::pretty_breaks(n = 20))
ggsave(output_pdf, width = 15,  height = 15,  units = c("in"), g6)