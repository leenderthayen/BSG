# - Define datasets known and used by BSG
#

# - ENSDF_099
bsg_add_dataset(
  NAME      ENSDF_099
  VERSION   201809
  URL       https://www.nndc.bnl.gov/ensarchivals/distributions/dist18/ensdf_180901_099.zip
  ENVVAR    ENSDF
  )

# - ENSDF_199
bsg_add_dataset(
  NAME      ENSDF_199
  VERSION   201809
  URL       https://www.nndc.bnl.gov/ensarchivals/distributions/dist18/ensdf_180901_199.zip
  ENVVAR    ENSDF
  )

# - ENSDF_299
bsg_add_dataset(
  NAME      ENSDF_299
  VERSION   201809
  URL       https://www.nndc.bnl.gov/ensarchivals/distributions/dist18/ensdf_180901_299.zip
  ENVVAR    ENSDF
  )

# - FRDM202
bsg_add_dataset(
  NAME      FRDM2012
  VERSION   2012
  URL       https://ars.els-cdn.com/content/image/1-s2.0-S0092640X1600005X-mmc1.zip
  ENVVAR    FRDM2012
  )

# - ChargeRadii
bsg_add_dataset(
  NAME      ChargeRadii
  VERSION   2016
  URL       https://journals.aps.org/prc/supplemental/10.1103/PhysRevC.94.064315/nuclear_charge_radii.txt
  ENVVAR    ChargeRadii
  )
