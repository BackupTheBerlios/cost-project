#===================================
# EJ, 24/09/2007
# lsData-class
#===================================

setClass("lsData",
	representation(
		desc="character",
		cc="data.frame"
	),
	prototype(
		desc="my stock",
		cc=data.frame(
			record_type="cc", 
			landing_country=NA, 
			vessel_flag=NA, 
			year=NA, 
			quarter=NA, 
			month=NA, 
			division=NA, 
			rectangle=NA, 
			species=NA, 
			landings_cat=NA, 
			comm_cat=NA, 
			fishing_act_nat=NA, 
			fishing_act_eu=NA, 
			unallocated_catch=NA,
			misreported_catch=NA,
			official_landings=NA,
			landings_multiplier=NA,
			currency=NA,
			landings_value=NA)		
	)
)
