#===================================
# EJ, 24/09/2007
# esData-class
#===================================

setClass("esData",
	representation(
		desc="character",
		es="data.frame"
	),
	prototype(
		desc="my stock",
		es=data.frame(
			record_type="es", 
			landing_country=NA, 
			vessel_flag=NA, 
			year=NA, 
			quarter=NA, 
			month=NA, 
			division=NA, 
			rectangle=NA, 
			fishing_act_nat=NA, 
			fishing_act_eu=NA, 
			n_trips=NA,
			n_fo=NA,
			effort_hours=NA,
			effort_kw_days=NA,
			effort_gt_days=NA,
			effort_days=NA)		
	)
)
