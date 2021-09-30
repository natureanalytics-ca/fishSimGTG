

#Life history object
setClass("LifeHistory",
         representation(
           title = "character",
           description = "character",
           speciesName = "character",
           L_units = "character",
           Walpha_units = "character",
           Linf = "numeric",
           K =  "numeric",
           L50 = "numeric",
           L95 = "numeric",
           M =  "numeric",
           MK = "numeric",
           LW_A = "numeric",
           LW_B = "numeric",
           Tmax = "numeric",
           Steep = "numeric")
)

#lbspr Object
setClass("LBSPRarray",
  representation(
    LifeHistory = "LifeHistory",
    sim = "list"
  )
)

