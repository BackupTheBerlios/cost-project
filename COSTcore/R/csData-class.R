#===================================
#
# EJ, 24/09/2007
# csData-class
#
#===================================

# 
# samp = sampling
# country = ctry
# flag =flg 
# vessel = vsl
# trip = trp
# number = num
# project = proj
# power = pwr
# fishing operation = fo
# method = meth
# depth = dep
# selectivity =sel
# device =dev
# commercial = comm
# weight = wt
# class = cls
# otolith = oto

setClass("csData",
	representation(
		desc="character",
		tr="data.frame",
		hh="data.frame",
		sl="data.frame",
		hl="data.frame",
		ca="data.frame"
	),
	prototype(
		desc="my stock",
		tr=data.frame(
			sampType=NA, # PK
			landCtry=NA, # PK
			vslFlgCtry=NA, # PK
			year=NA, # PK
			proj=NA, # PK
			trpNum=NA, # PK
			vslLen=NA, 
			vslPwr=NA, 
			vslSize=NA, 
			vsType=NA, 
			foNum=NA, 
			daysAtSea=NA, 
			vslId=NA, 
			sampCtry=NA, 
			sampMeth=NA),
		hh=data.frame(
			sampType=NA, # FK
			landCtry=NA, # FK
			vslFlgCtry=NA, # FK
			year=NA, # FK
			proj=NA, # FK
			trpNum=NA, # FK
			staNum=NA, # PK
			foVal=NA,
			aggLev=NA,
			date=NA,
			timeShot=NA,
			foDur=NA,
			latIni=NA,
			lonIni=NA,
			latFin=NA,
			lonFin=NA,
			area=NA,
			rect=NA,
			foDep=NA,
			waterDep=NA,
			foCatNat=NA,
			foCatEu5=NA,
			foCatEu6=NA,
			gear=NA,
			meshSize=NA,
			selDev=NA,
			meshSizeSelDev=NA),
		sl=data.frame(
			sampType=NA, # FK
			landCtry=NA, # FK
			vslFlgCtry=NA, # FK
			year=NA, # FK
			proj=NA, # FK
			trpNum=NA, # FK
			staNum=NA, # FK
			spp=NA, # PK 
			sex=NA, # PK
			catchCat=NA, # PK 
			landCat=NA, # PK 
			commCatScl=NA, # PK
			commCat=NA, # PK
			subSampCat=NA, # PK
			valCode=NA, 
			wt=NA, 
			subSampWt=NA, 
			lenCode=NA),
		hl=data.frame(
			sampType=NA, # FK
			landCtry=NA, # FK
			vslFlgCtry=NA, # FK
			year=NA, # FK
			proj=NA, # FK
			trpNum=NA, # FK
			staNum=NA, # FK
			spp=NA, # PK 
			sex=NA, # PK
			catchCat=NA, # PK 
			landCat=NA, # PK 
			commCatScl=NA, # PK
			commCat=NA, # PK
			subSampCat=NA, # PK
			lenCls=NA, 
			lenNum=NA),
		ca=data.frame(
			sampType=NA, # FK
			landCtry=NA, # FK
			vslFlgCtry=NA, # FK
			year=NA, # FK
			proj=NA, # FK
			trpNum=NA, # FK
			staNum=NA, # PK (optional)
			spp=NA, # PK 
			sex=NA, # PK
			catchCat=NA, # PK 
			landCat=NA, # PK 
			commCatScl=NA, # PK
			commCat=NA, # PK
			stock=NA, # PK
			area=NA, # PK
			rect=NA, # PK
			lenCode=NA, # PK
			lenCls=NA, # PK
			age=NA, # PK
			plusGrp=NA, # PK
			otoWt=NA,
			otoSide=NA,
			indWt=NA,
			matScale=NA,
			matStage=NA,
			num=NA,
			fishId=NA) # PK 
	)
)
