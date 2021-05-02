#this installs TCGAbiolinks if needed
if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)



install.packages("survival")
install.packages("survminer")
install.packages("arsenal")

#if you get an error saying you need a CRAN mirror, it is with the above. Opening an R terminal
# and manually installing the packages worked for me

 library(survival)
 library(survminer)
 library(arsenal)   

 
clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", file.type="xml")
GDCdownload( clin_query ) #should only need this command once. This downloads the files onto your system.
clinic <- GDCprepare_clinic(clin_query, clinical.info="patient")  #these download clinical data into clinic matrix
names(clinic)[names(clinic) == "days_to_last_followup"] = "days_to_last_follow_up"   # fixes missing
# underscore, would interupt TCGA_analyze looking for that column otherwise


age_clinical = clinic$age_at_initial_pathologic_diagnosis
clinic$age_category = ifelse(age_clinical <40, "Young", ifelse(age_clinical>=60, "Old", "Mid"))
#adds column classifying each patient as young, mid or old.

subtypes <- TCGAquery_subtype(tumor = "BRCA")
#grabs the database with PAM50 subtypes, clinic doesn't actually have

barcodes_subtypes = subtypes$patient
PAM50 = rep("NA",length(clinic$bcr_patient_barcode))
barcodes_clinic = clinic$bcr_patient_barcode
counter2=0
for (barcode1 in barcodes_clinic){
	counter = 0
	counter2= counter2 +1
	for (barcode2 in barcodes_subtypes){
		counter= counter +1
		if( barcode1==barcode2){
			PAM50[counter2] = subtypes$BRCA_Subtype_PAM50[counter]
		}
		
	}
}

clinic$PAM50 = PAM50

#this grabs the subtype column and puts it in clinic. Doing TCGA_analyze from subtypes gives errors that 
#I didn't want to try deciphering
#the nested for loops are needed cause subtypes is smaller. Doesn't have all the same patients, so need 
#verify barcodes


NA_list = c()

young_list = c()
counter = 0
for (tumor in PAM50){
	counter= counter +1
	if (tumor =="NA"){
		NA_list =c(NA_list,counter)
	}
}
#NA list now holds row numbers of Patients with no listed subtype
counter = 0
age_groups = clinic$age_category
for (age in age_groups){
        counter= counter +1
        if (age =="Young"){
                young_list =c(young_list,counter)
        }
}
#young list now holds row numbers of patients that are young

copy_clinic = clinic[-NA_list,] #copy is clinic but NA subtype rows excluded
young_clinic = copy_clinic[young_list,] #young is copy_clinic but just young (also no NA) 
old_clinic = copy_clinic[-young_list,] #old is copy_clinic but excluding young (also no NA)

TCGAanalyze_survival( copy_clinic, "PAM50", filename="figures/PAM50_survival.pdf")

TCGAanalyze_survival( young_clinic, "PAM50", filename="figures/young_PAM50_survival.pdf")

TCGAanalyze_survival( old_clinic, "PAM50", filename="figures/older_PAM50_survival.pdf")

#TCGA_analyze grabs several columns and provided column to make Kaplan Meier survival plots



jpeg("figures/distributetest.jpg")
barplot(table(copy_clinic$PAM50), main = "Distribution of Subtypes for Patients in TCGA Database",
	xlab = "PAM50 Subtype", ylab ="Number of Tumors")
dev.off()

#barplot makes a bargraph. Here shows distribution of subtypes. Using modified clinics to focus on young or old


jpeg("figures/youngdistributetest.jpg")
barplot(table(young_clinic$PAM50), main = "Subtype distribution for young patients(<40 years)",
        xlab = "PAM50 Subtype", ylab ="Number of Tumors")
dev.off()

jpeg("figures/olderdistributetest.jpg")
barplot(table(old_clinic$PAM50), main = "Subtype distribution for older patients(>=40 years)",
        xlab = "PAM50 Subtype", ylab ="Number of Tumors")
dev.off()
