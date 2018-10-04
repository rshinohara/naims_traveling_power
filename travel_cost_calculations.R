############################################################################################################
# Travel Cost Calculations for NAIMS Traveling Study
# Russell Shinohara
# October 2018 Submission
############################################################################################################

###
# MISCELLANEOUS EXPENSES
###
# https://public.csr.nih.gov/ReviewerResources/TravelAndExpenses/Documents/NonFederalPeerReviewTravelGuidelines.pdf
# Accessed December 19, 2017 - SREA Reimbursement - $235 for miscellaneous expenses (taxis, shuttles, luggage fees, 
# internet Wi-Fi, airport parking fees, tolls, mileage reimbursement to and from the home destination, etc.)
misc.cost<-235

###
# TRANSIT COSTS
###
sites<-c('PEN','CHP','COR','NIH','HOP','VAN','BUF','BWH','MGH','UVM')
cities<-c('PHL','PHL','NYC','WAS','BAL','NAS','BUF','BOS','BOS','BUR')

cities.unique<-c('PHL','NYC','WAS','BAL','NAS','BUF','BOS','BUR')
travel.costs.cities<-read.csv('~/Dropbox/Documents/Grants/NIH/NAIMS_Traveling_R01_2017/third_submission/travel_costs.csv',header=FALSE)
travel.costs.cities<-t(matrix(as.numeric(levels(unlist(travel.costs.cities))[unlist(travel.costs.cities)]),nrow=8))
rownames(travel.costs.cities)<-cities.unique
colnames(travel.costs.cities)<-cities.unique

###
# HOTEL COSTS
###

hotel.costs<-c(188.33,199.00,189.00,128.67,227.33,147.00,156.33,252.33)

############################################################################################################
############################################################################################################
# PHASE I COSTS
############################################################################################################
############################################################################################################
n.1<-3 # number of subjects in Phase I

# Participants will be recruited at random from the sites, and then sent to all sites in three steps:
phase1.travelcost<-rep(NA,length(hotel.costs))

# For a Philadelphia-based participant
hotel.costs.this<-2*sum(hotel.costs[-1])
travel.costs<-rep(NA,3)
# Travel PHL->BUF->BUR->PHL
travel.costs[1]<-travel.costs.cities[1,6]+travel.costs.cities[6,8]+travel.costs.cities[8,1]
# Travel PHL->BOS->NYC->PHL
travel.costs[2]<-travel.costs.cities[1,7]+travel.costs.cities[7,2]+travel.costs.cities[2,1]
# Travel PHL->BAL->DC->NAS->PHL
travel.costs[3]<-travel.costs.cities[1,4]+travel.costs.cities[4,3]+travel.costs.cities[3,5]+travel.costs.cities[5,1]
phase1.travelcost[1]<-hotel.costs.this+sum(travel.costs)+misc.cost*7

# For a New York-based participant
hotel.costs.this<-2*sum(hotel.costs[-2])
travel.costs<-rep(NA,3)
# Travel NYC->BUF->BUR->NYC
travel.costs[1]<-travel.costs.cities[2,6]+travel.costs.cities[6,8]+travel.costs.cities[8,2]
# Travel NYC->BOS->PHL->NYC
travel.costs[2]<-travel.costs.cities[2,7]+travel.costs.cities[7,1]+travel.costs.cities[1,2]
# Travel NYC->BAL->DC->NAS->NYC
travel.costs[3]<-travel.costs.cities[2,4]+travel.costs.cities[4,3]+travel.costs.cities[3,5]+travel.costs.cities[5,2]
phase1.travelcost[2]<-hotel.costs.this+sum(travel.costs)+misc.cost*7

# For a DC-based participant
hotel.costs.this<-2*sum(hotel.costs[-3])
travel.costs<-rep(NA,3)
# Travel DC->BUF->BUR->DC
travel.costs[1]<-travel.costs.cities[3,6]+travel.costs.cities[6,8]+travel.costs.cities[8,3]
# Travel DC->BOS->NYC->DC
travel.costs[2]<-travel.costs.cities[3,7]+travel.costs.cities[7,2]+travel.costs.cities[2,3]
# Travel DC->PHL->BAL->NAS->DC
travel.costs[3]<-travel.costs.cities[3,1]+travel.costs.cities[1,4]+travel.costs.cities[4,5]+travel.costs.cities[5,3]
phase1.travelcost[3]<-hotel.costs.this+sum(travel.costs)+misc.cost*7

# For a Baltimore-based participant
hotel.costs.this<-2*sum(hotel.costs[-4])
travel.costs<-rep(NA,3)
# Travel BAL->BUF->BUR->BAL
travel.costs[1]<-travel.costs.cities[4,6]+travel.costs.cities[6,8]+travel.costs.cities[8,4]
# Travel BAL->BOS->NYC->BAL
travel.costs[2]<-travel.costs.cities[4,7]+travel.costs.cities[7,2]+travel.costs.cities[2,4]
# Travel BAL->PHL->DC->NAS->BAL
travel.costs[3]<-travel.costs.cities[4,1]+travel.costs.cities[1,3]+travel.costs.cities[3,5]+travel.costs.cities[5,4]
phase1.travelcost[4]<-hotel.costs.this+sum(travel.costs)+misc.cost*7

# For a Nashville-based participant
hotel.costs.this<-2*sum(hotel.costs[-5])
travel.costs<-rep(NA,3)
# Travel NAS->BUF->BUR->NAS
travel.costs[1]<-travel.costs.cities[5,6]+travel.costs.cities[6,8]+travel.costs.cities[8,5]
# Travel NAS->BOS->NYC->NAS
travel.costs[2]<-travel.costs.cities[5,7]+travel.costs.cities[7,2]+travel.costs.cities[2,5]
# Travel NAS->PHL->DC->BAL->NAS
travel.costs[3]<-travel.costs.cities[5,1]+travel.costs.cities[1,3]+travel.costs.cities[3,4]+travel.costs.cities[4,5]
phase1.travelcost[5]<-hotel.costs.this+sum(travel.costs)+misc.cost*7

# For a Buffalo-based participant
hotel.costs.this<-2*sum(hotel.costs[-6])
travel.costs<-rep(NA,3)
# Travel BUF->BUR->BOS->BUF
travel.costs[1]<-travel.costs.cities[6,8]+travel.costs.cities[8,7]+travel.costs.cities[7,6]
# Travel BUF->NYC->PHL->BUF
travel.costs[2]<-travel.costs.cities[6,2]+travel.costs.cities[2,1]+travel.costs.cities[1,6]
# Travel BUF->BAL->DC->NAS->BUF
travel.costs[3]<-travel.costs.cities[6,4]+travel.costs.cities[4,3]+travel.costs.cities[3,5]+travel.costs.cities[5,7]
phase1.travelcost[6]<-hotel.costs.this+sum(travel.costs)+misc.cost*7

# For a Boston-based participant
hotel.costs.this<-2*sum(hotel.costs[-7])
travel.costs<-rep(NA,3)
# Travel BOS->BUF->BUR->BOS
travel.costs[1]<-travel.costs.cities[7,6]+travel.costs.cities[6,8]+travel.costs.cities[8,7]
# Travel BOS->NYC->PHL->BOS
travel.costs[2]<-travel.costs.cities[7,2]+travel.costs.cities[2,1]+travel.costs.cities[1,7]
# Travel BOS->BAL->DC->NAS->BOS
travel.costs[3]<-travel.costs.cities[7,4]+travel.costs.cities[4,3]+travel.costs.cities[3,5]+travel.costs.cities[5,7]
phase1.travelcost[7]<-hotel.costs.this+sum(travel.costs)+misc.cost*7

# For a Burlington-based participant
hotel.costs.this<-2*sum(hotel.costs[-8])
travel.costs<-rep(NA,3)
# Travel BUR->BOS->BUF->BUR
travel.costs[1]<-travel.costs.cities[8,7]+travel.costs.cities[7,6]+travel.costs.cities[6,8]
# Travel BUR->NYC->PHL->BUR
travel.costs[2]<-travel.costs.cities[8,2]+travel.costs.cities[2,1]+travel.costs.cities[1,8]
# Travel BUR->BAL->DC->NAS->BUR
travel.costs[3]<-travel.costs.cities[8,4]+travel.costs.cities[4,3]+travel.costs.cities[3,5]+travel.costs.cities[5,8]
phase1.travelcost[8]<-hotel.costs.this+sum(travel.costs)+misc.cost*7

#These are travel cost estimates for one time point only
mean(phase1.travelcost)

## PHASE I TRAVEL COSTS = ~ 5400 per participant per timepoint
mean(phase1.travelcost)*2*n.1
# = 32,200 for two time points


############################################################################################################
############################################################################################################
# PHASE II COSTS
############################################################################################################
############################################################################################################

n.2<-30 # number of subjects in Phase II
n.scans<-3 # number of imaging visits for each subject
n.subj.persite<-n.2/length(sites) # number of imaging visits for each subject

sample.sites<-function(n,origin.site=NULL) {
  possible.sites<-1:length(sites)
  if (!is.null(origin.site)) possible.sites<-possible.sites[-which(cities==cities[origin.site])]
	j<-sample(possible.sites,n,replace=FALSE)
	list(sites=sites[j],cities=cities[j],site.numbers=j)
}

#includes ground transportation/meals
get.travel.costs<-function(dests,origin) {
  total.travel.costs<-travel.costs.cities[origin,dests[1]]+travel.costs.cities[dests[length(dests)],origin]
  if (length(dests)>1) for (l in 2:(length(dests))) total.travel.costs<-total.travel.costs+travel.costs.cities[dests[l-1],dests[l]]
  total.travel.costs+length(dests)*misc.cost+2*sum(hotel.costs[dests])
}

#These are travel cost estimates for one time point only
estimate.total.travel.costs<-function() {
  total.travel.costs<-0
  for (l in 1:length(sites)) {
    for (i in 1:n.scans) {
      sampled<-sample.sites(n.scans-1,origin.site=l)
      unique.site.numbers<-sampled$site.numbers[!duplicated(sampled$cities)]
      dests<-which(rownames(travel.costs.cities)%in%sampled$cities[!duplicated(sampled$cities)]) 
      total.travel.costs<-total.travel.costs+get.travel.costs(dests,origin=which(rownames(travel.costs.cities)==cities[l]))
    }
  }
  total.travel.costs
}

## PHASE II TRAVEL COSTS = ~ 46,000 per timepoint
set.seed(32542)
mean(replicate(1000,estimate.total.travel.costs()))*2
# = 92,000 for two time points

############################################################################################################
############################################################################################################
# TOTAL COST SUMMARY
############################################################################################################
############################################################################################################

# TOTAL TRAVEL COSTS
set.seed(32542); mean(phase1.travelcost)*2*n.1+mean(replicate(1000,estimate.total.travel.costs()))*2
# = 124,100

