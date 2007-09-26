#===================================
#
# EJ, 24/09/2007
# ceData-class
#
#===================================

setClass("ceData",
	representation(
		desc="character",
		ce="data.frame"
	),
	prototype(
		desc="my stock",
		ce=data.frame(
			vslFlgCtry=NA, # PK
			year=NA, # PK
			quarter=NA, # PK 
			month=NA, # PK
			area=NA, # PK
			rectangle=NA, # PK 
			landCat=NA, # PK 
			foCatNat=NA, # PK
			foCatEu5=NA, # PK
			foCatEu6=NA, # PK
			trpNum=NA,
			foNum=NA,
			foDur=NA,
			effKwDays=NA,
			effGtDays=NA,
			daysAtSea=NA)		
	)
)
