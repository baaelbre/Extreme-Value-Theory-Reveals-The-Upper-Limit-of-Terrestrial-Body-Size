# Install and load necessary packages
if (!require("paleobioDB")) install.packages("paleobioDB")
library(paleobioDB)

# Search for sauropods
sauropods <- pbdb_occurrences(base_name = "Sauropoda", show = c("phylo", "coords", "ident"))
# Get occurrence numbers
occurrence_numbers <- sauropods$oid

all_measurements <- pbdb_measurements(occ_id = occurrence_numbers, show='spec')

# TODO: from here onwards ...
# Filter for femur measurements
femur_circ <- all_measurements[which(all_measurements$smp == "femur" & all_measurements$mty == 'circumference'), ]
femur_length <- all_measurements[which(all_measurements$smp == "femur" & all_measurements$mty == 'length'), ]

# Filter for humerus measurements
humerus_circ <- all_measurements[which(all_measurements$smp == "humerus" & all_measurements$mty == 'circumference'), ]
humerus_length <- all_measurements[which(all_measurements$smp == "humerus" & all_measurements$mty == 'length'), ]
