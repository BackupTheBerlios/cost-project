#===================================
#
# EJ, 24/09/2007
# clData-class
#
#===================================

setClass("clData",
	representation(
		desc="character",
		cl="data.frame"
	),
	prototype(
		desc="my stock",
		cl=data.frame(
			landCtry=NA, # PK
			vslFlgCtry=NA, # PK
			year=NA, # PK
			quarter=NA, # PK 
			month=NA, # PK
			area=NA, # PK
			rectangle=NA, # PK 
			spp=NA, # PK 
			landCat=NA, # PK 
			commCatScl=NA, # PK
			commCat=NA, # PK
			foCatNat=NA, # PK
			foCatEu5=NA, # PK
			foCatEu6=NA, # PK
			unallocCatchWt=NA,
			misRepCatchWt=NA,
			landWt=NA,
			landMult=NA,
			landValue=NA)		
	)
)
