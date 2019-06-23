library(tidyverse)

read_obj <- function(inobj, flag_fsep=0) {
  cat("read from",inobj, "\n")
  tmp <- read_csv(inobj, col_names = FALSE, comment = "#", col_types=cols(X1=col_character()))
  
  tmpv <- tmp %>%
    dplyr::filter(str_detect(X1,"^v ")) %>%
    separate("X1", c('K', 'X','Y','Z'), sep=" ") %>%
    select(-K)
  tmpv$X <- parse_double(tmpv$X)
  tmpv$Y <- parse_double(tmpv$Y)
  tmpv$Z <- parse_double(tmpv$Z)
  vert <<- tmpv
  #write_csv(tmpv, sub(".obj","_obj_vert.txt",inobj))
  
  tmpf <- tmp %>%
    dplyr::filter(str_detect(X1,"^f ")) %>%
    separate("X1", c('K', 'F1','F2','F3'), sep=" ") %>%
    select(-K)
  # --- extract v from face data ---
  if(flag_fsep == 0) {
    tmpf$F1 = parse_integer(sub("/.*","",tmpf$F1))
    tmpf$F2 = parse_integer(sub("/.*","",tmpf$F2))
    tmpf$F3 = parse_integer(sub("/.*","",tmpf$F3))
    face <<- tmpf
    #write_csv(tmpf, sub(".obj","_obj_face.txt",inobj))
  } else {
    names(tmpf) <- c("F1a", "F2a", "F3a")
    tmpff <- tmpf %>%
      separate("F1a", c("F1","F1vt","F1vn"), sep="/",fill='right') %>%
      separate("F2a", c("F2","F2vt","F2vn"), sep="/",fill='right') %>%
      separate("F3a", c("F3","F3vt","F3vn"), sep="/",fill='right')
    tmpff$F1   = parse_integer(tmpff$F1)
    tmpff$F1vt = parse_integer(tmpff$F1vt)
    tmpff$F1vn = parse_integer(tmpff$F1vn)
    tmpff$F2   = parse_integer(tmpff$F2)
    tmpff$F2vt = parse_integer(tmpff$F2vt)
    tmpff$F2vn = parse_integer(tmpff$F2vn)
    tmpff$F3   = parse_integer(tmpff$F3)
    tmpff$F3vt = parse_integer(tmpff$F3vt)
    tmpff$F3vn = parse_integer(tmpff$F3vn)
    face <<- tmpff
  }
}

# === Normal Vector ===
# pos1,pos2,pos3 : vetices of face
NormalVector3D <- function(pos1, pos2, pos3) {
  vec1 <- SubVector3D( pos2, pos1 )
  vec2 <- SubVector3D( pos3, pos1 )
  return(CrossProduct3D( vec1, vec2 ))
}

# === Subtraction of 3D Vector ===
SubVector3D <- function( vec1, vec2 ) {
  return(c(vec1[1] - vec2[1],
           vec1[2] - vec2[2],
           vec1[3] - vec2[3]))
}

# === inner product ===
DotProduct3D <- function( vec1, vec2 ) {
  return ( vec1[1] * vec2[1] + vec1[2] * vec2[2] + vec1[3] * vec2[3] )
}

#  === cross product ===
CrossProduct3D <- function(vec1, vec2) {
  return(c(vec1[2] * vec2[3] - vec1[3] * vec2[2],
           vec1[3] * vec2[1] - vec1[1] * vec2[3],
           vec1[1] * vec2[2] - vec1[2] * vec2[1]))
}

# ===  cross point between line and plane ===
# pos1 : a point(vertex) on line
# vec1 : a direction vector of line
# pos2 : a point(vertex) on plane
# vec2 : a normal vector of plane
LinePolyCrossPoint <- function(pos1s, pos1e, pos2, vec2) {
  vec1 <- SubVector3D( pos1e, pos1s)
  #//PrintPoint3(pos1s)
  #//PrintPoint3(vec1)
  
  vec <- SubVector3D( pos2, pos1s)
  a <- DotProduct3D( vec, vec2 )
  b <- DotProduct3D( vec1, vec2 )
  
  pos <- c(0,0,0)
  # --- Since they are parallel, there is no intersection ---
  if( b == 0.0 ) {
    return(list(0,pos))
  }
  # --- Find the parameter t ---
  t <- a / b
  #//fprintf(stderr, "t = %f\n", t)
  
  # --- Since it is out of range, there is no intersection ---
  if( t < 0.0 || 1.0 < t ) {
    return(list(2,pos))
  }
  # --- Substitute into the equation of the straight line and find the intersection point ---
  pos <- c(pos1s[1] + vec1[1] * t,
           pos1s[2] + vec1[2] * t,
           pos1s[3] + vec1[3] * t)
  
  # --- With intersection ---
  return(list(1,pos))
}

PrintPoint3 <- function( pos ) {
  sprintf("check : %f %f %f\n", pos[1], pos[2], pos[3])
}

# === make profile ===
profile <- function(pid, prof_pos, prof_nvec, out_prof) {
  #prof_pos=c(0,0,0)
  #prof_nvec=c(0,1,0)
  
  if(pid == 1) {
    # prof_file <- file(out_prof, "w")
    cat("profile",pid," ","create",out_prof," ")
  } else {
    # prof_file <- file(out_prof, "a")
    cat("profile",pid," ","append",out_prof," ")
  }
  
  posSave <- list(c(0,0),c(0,0),c(0,0))
  names(posSave) <- c("X","Y","Z")
  
  flag_prlist_first <- 1
  prlist <- data.frame(ID=1, X=0, Y=0, Z=0)
  prlist <- prlist[-1,]
  
  # --- convert to vector for speed up ---
  f1 <- face$F1
  f2 <- face$F2
  f3 <- face$F3
  vx <- vert$X
  vy <- vert$Y
  vz <- vert$Z
  
  # --- main loop ---
  cat("|")
  for(i in 1:length(f1)) {
    if((i %% 100000) == 0) {
      cat("=")
    }
    
    #t1 <- face$F1[i]
    t1 <- f1[i]
    t2 <- f2[i]
    t3 <- f3[i]
    
    x <- c(vx[t1],vx[t2],vx[t3])
    y <- c(vy[t1],vy[t2],vy[t3])
    z <- c(vz[t1],vz[t2],vz[t3])
    
    # --- pre-selection ---
    if((prof_nvec[1] == 0) && (prof_nvec[2] == 1) && (prof_nvec[3] == 0)) {
      if(prof_pos[2] < (min(y)-0.1) || (max(y)+0.1) < prof_pos[2] ) {
        next()
      }
    } else if((prof_nvec[1] == 1) && (prof_nvec[2] == 0) && (prof_nvec[3] == 0)) {
      if(prof_pos[1] < (min(x)-0.1) || (max(x)+0.1) < prof_pos[1] ) {
        next()
      }
    } else if((prof_nvec[1] == 0) && (prof_nvec[2] == 0) && (prof_nvec[3] == 1)) {
      if(prof_pos[3] < (min(z)-0.1) || (max(z)+0.1) < prof_pos[3] ) {
        next()
      }
    }
    
    cross_num <- 0
    for(ts in 1:3) {
      # --- start point ---
      pos1s <- c( x[ts], y[ts], z[ts] )
      
      # --- end point ---
      te <- ts + 1
      if(te == 4) {
        te <- 1
      }
      pos1e <- c( x[te], y[te], z[te] )
      
      # --- cross point between line and plane ---
      ret <- LinePolyCrossPoint(pos1s, pos1e, prof_pos, prof_nvec)
      flag <- ret[[1]]
      pos <- ret[[2]]
      if(flag == 1) {
        cross_num <- cross_num + 1
        posSave$X[cross_num] <- round(pos[1],digit=6) #format(pos[1], nsmall = 6)
        posSave$Y[cross_num] <- round(pos[2],digit=6) #format(pos[2], nsmall = 6)
        posSave$Z[cross_num] <- round(pos[3],digit=6) #format(pos[3], nsmall = 6)
      }
      
      if(cross_num == 2) {
        #write string is return
        #write(c(pid, posSave$X[1],posSave$Y[1],posSave$Z[1]), file=prof_file, sep=" ")
        #write(c(pid, posSave$X[2],posSave$Y[2],posSave$Z[2]), file=prof_file, sep=" ")
        prlist1 <- data.frame(ID=c(pid,pid), X=posSave$X[1:2], Y=posSave$Y[1:2], Z=posSave$Z[1:2])
        prlist <- rbind(prlist, prlist1)
      }
    }
  }
  gprlist <<- prlist
  if(pid == 1) {
    write_csv(prlist, out_prof)
  } else {
    write_csv(prlist, out_prof, append = TRUE)
  }
  
  cat("|\n")
  #close(prof_file)
}

# === bouding box ===
bbox <- function(out_bbox) {
  cat("crete",out_bbox,"\n")
  bbox_file <- file(out_bbox, "w")
  # --- convert to vector for speed up ---
  vx <- vert$X
  vy <- vert$Y
  vz <- vert$Z
  xmin <- min(vx)
  xmax <- max(vx)
  ymin <- min(vy)
  ymax <- max(vy)
  zmin <- min(vz)
  zmax <- max(vz)
  # write("#Box dimensions", file=bbox_file, sep=" ")
  # write(sprintf("X: %f",xmax-xmin), file=bbox_file, sep=" ")
  # write(sprintf("Y: %f",ymax-ymin), file=bbox_file, sep=" ")
  # write(sprintf("Z: %f",zmax-zmin), file=bbox_file, sep=" ")
  # 
  # write("", file=bbox_file, sep=" ")
  # write("#Box center", file=bbox_file, sep=" ")
  # write(sprintf("X: %f",(xmax+xmin)/2), file=bbox_file, sep=" ")
  # write(sprintf("Y: %f",(ymax+ymin)/2), file=bbox_file, sep=" ")
  # write(sprintf("Z: %f",(zmax+zmin)/2), file=bbox_file, sep=" ")
  # 
  # write("", file=bbox_file, sep=" ")
  # write("#Bounding Box", file=bbox_file, sep=" ")
  # write(sprintf("xmin: %f",xmin), file=bbox_file, sep=" ")
  # write(sprintf("xmax: %f",xmax), file=bbox_file, sep=" ")
  # write(sprintf("ymin: %f",ymin), file=bbox_file, sep=" ")
  # write(sprintf("ymax: %f",ymax), file=bbox_file, sep=" ")
  # write(sprintf("zmin: %f",zmin), file=bbox_file, sep=" ")
  # write(sprintf("zmax: %f",zmax), file=bbox_file, sep=" ")
  
  write(c(xmin, ymin, zmin), file=bbox_file, sep=" ")
  write(c(xmax, ymin, zmin), file=bbox_file, sep=" ")
  write(c(xmax, ymax, zmin), file=bbox_file, sep=" ")
  write(c(xmin, ymax, zmin), file=bbox_file, sep=" ")
  write(c(xmin, ymin, zmax), file=bbox_file, sep=" ")
  write(c(xmax, ymin, zmax), file=bbox_file, sep=" ")
  write(c(xmax, ymax, zmax), file=bbox_file, sep=" ")
  write(c(xmin, ymax, zmax), file=bbox_file, sep=" ") 
  close(bbox_file)
}

# === multiple profiles from y ===
profiles_y <- function(ndiv, prof_file) {
  prof_nvec <- c(0,1,0)
  # --- convert to vector for speed up ---
  vx <- vert$X
  vy <- vert$Y
  vz <- vert$Z
  # --- maxmin ---
  ymin <- min(vy)
  ymax <- max(vy)
  # --- loop ---
  ystep <- (ymax - ymin) / ndiv
  for(pid in c(1:(ndiv-1))) {
    #profile(c(0,ymin+pid*ystep,0), prof_nvec, sprintf("%s_%03d.txt", prof_base, pid))
    profile(pid, c(0,ymin+pid*ystep,0), prof_nvec, prof_file)
  }
}
profiles_x <- function(ndiv, prof_file) {
  prof_nvec <- c(1,0,0)
  # --- convert to vector for speed up ---
  vx <- vert$X
  vy <- vert$Y
  vz <- vert$Z
  # --- maxmin ---
  xmin <- min(vx)
  xmax <- max(vx)
  # --- loop ---
  xstep <- (xmax - xmin) / ndiv
  for(pid in c(1:(ndiv-1))) {
    profile(pid, c(xmin+pid*xstep,0,0), prof_nvec, prof_file)
  }
}
profiles_z <- function(ndiv, prof_file) {
  prof_nvec <- c(0,0,1)
  # --- convert to vector for speed up ---
  vx <- vert$X
  vy <- vert$Y
  vz <- vert$Z
  # --- maxmin ---
  zmin <- min(vz)
  zmax <- max(vz)
  # --- loop ---
  zstep <- (zmax - zmin) / ndiv
  for(pid in c(1:(ndiv-1))) {
    profile(pid, c(0,0,zmin+pid*zstep), prof_nvec, prof_file)
  }
}

read_objs <- function(dir1) {
  objs <- list.files(dir1, pattern=".obj", full.names=T)
  for(fpath in objs) {
    # === read obj ===
    # readobj(fpath)
    read_obj(fpath)
  }
}

objs_profies_xyz <- function(dir1, ndiv) {
  objs <- list.files(dir1, pattern=".obj$", full.names=T)
  for(fpath in objs) {
    # === read obj ===
    # readobj(fpath)
    read_obj(fpath)
    
    # === bbox ===
    bbox_file <- sub(".obj","_bbox.txt",fpath)
    bbox(bbox_file)
    
    # === profile ===
    profiles_y(ndiv,sub(".obj","_prof_y.txt",fpath))
    profiles_x(ndiv,sub(".obj","_prof_x.txt",fpath))
    profiles_z(ndiv,sub(".obj","_prof_z.txt",fpath))
  }
}

read_profiles_csv <- function(prof_file) {
  prof <- read_csv(prof_file,
                   col_types = cols(
                     ID = col_integer(),
                     X = col_double(),
                     Y = col_double(),
                     Z = col_double()
                   )
  )
  prof
}
read_profiles <- function(prof_file) {
  prof <<- read_table2(prof_file,
                       col_names = c("ID", "X", "Y", "Z"),
                       col_types = cols(
                         ID = col_integer(),
                         X = col_double(),
                         Y = col_double(),
                         Z = col_double())
  )
}

profiles_y_bbox <- function(prof, prof_bbox_file, flag_1line) {
  cat("write",prof_bbox_file,"\n")
  id_max <- max(prof$ID)
  for (i in c(1:id_max)) {
    p1 <- prof %>% dplyr::filter(ID==i)
    xx1 <- p1 %>% dplyr::filter(X==min(p1$X))
    xx2 <- p1 %>% dplyr::filter(X==max(p1$X))
    zz1 <- p1 %>% dplyr::filter(Z==min(p1$Z))
    zz2 <- p1 %>% dplyr::filter(Z==max(p1$Z))
    b <- xx1[1,] #becase result is multi, select first data.
    if(flag_1line == 1) {
      b <- cbind(b, xx2[1,-1])
      b <- cbind(b, zz1[1,-1])
      b <- cbind(b, zz2[1,-1])
      bbox <- b
      names(bbox) <- c("ID","x1","y1","z1","x2","y2","z2","x3","y3","z3","x4","y4","z4")
    } else {
      b <- rbind(b, xx2[1,])
      b <- rbind(b, zz1[1,])
      b <- rbind(b, zz2[1,])
      bb <- c("xx1","xx2","zz1","zz2")
      bbox <- data.frame(ID=b$ID, BBOX=bb, X=b$X, Y=b$Y, Z=b$Z)
    }
    if(i == 1) {
      write_csv(bbox,path = prof_bbox_file)
    } else {
      write_csv(bbox,path = prof_bbox_file, append = TRUE)
    }
  }
}

profiles_x_bbox <- function(prof, prof_bbox_file, flag_1line) {
  cat("write",prof_bbox_file,"\n")
  id_max <- max(prof$ID)
  for (i in c(1:id_max)) {
    p1 <- prof %>% dplyr::filter(ID==i)
    yy1 <- p1 %>% dplyr::filter(Y==min(p1$Y))
    yy2 <- p1 %>% dplyr::filter(Y==max(p1$Y))
    zz1 <- p1 %>% dplyr::filter(Z==min(p1$Z))
    zz2 <- p1 %>% dplyr::filter(Z==max(p1$Z))
    b <- yy1[1,] #becase result is multi, select first data.
    if(flag_1line == 1) {
      b <- cbind(b, yy2[1,-1])
      b <- cbind(b, zz1[1,-1])
      b <- cbind(b, zz2[1,-1])
      bbox <- b
      names(bbox) <- c("ID","x1","y1","z1","x2","y2","z2","x3","y3","z3","x4","y4","z4")
    } else {
      b <- rbind(b, yy2[1,])
      b <- rbind(b, zz1[1,])
      b <- rbind(b, zz2[1,])
      bb <- c("yy1","yy2","zz1","zz2")
      bbox <- data.frame(ID=b$ID, BBOX=bb, X=b$X, Y=b$Y, Z=b$Z)
    }
    if(i == 1) {
      write_csv(bbox,path = prof_bbox_file)
    } else {
      write_csv(bbox,path = prof_bbox_file, append = TRUE)
    }
  }
}

profiles_z_bbox <- function(prof, prof_bbox_file, flag_1line) {
  cat("write",prof_bbox_file,"\n")
  id_max <- max(prof$ID)
  for (i in c(1:id_max)) {
    p1 <- prof %>% dplyr::filter(ID==i)
    xx1 <- p1 %>% dplyr::filter(X==min(p1$X))
    xx2 <- p1 %>% dplyr::filter(X==max(p1$X))
    yy1 <- p1 %>% dplyr::filter(Y==min(p1$Y))
    yy2 <- p1 %>% dplyr::filter(Y==max(p1$Y))
    b <- xx1[1,] #becase result is multi, select first data.
    if(flag_1line == 1) {
      b <- cbind(b, xx2[1,-1])
      b <- cbind(b, yy1[1,-1])
      b <- cbind(b, yy2[1,-1])
      bbox <- b
      names(bbox) <- c("ID","x1","y1","z1","x2","y2","z2","x3","y3","z3","x4","y4","z4")
    } else {
      b <- rbind(b, xx2[1,])
      b <- rbind(b, yy1[1,])
      b <- rbind(b, yy2[1,])
      bb <- c("xx1","xx2","yy1","yy2")
      bbox <- data.frame(ID=b$ID, BBOX=bb, X=b$X, Y=b$Y, Z=b$Z)
    }
    if(i == 1) {
      write_csv(bbox,path = prof_bbox_file)
    } else {
      write_csv(bbox,path = prof_bbox_file, append = TRUE)
    }
  }
}

profies_xyz_bbox <- function(srcDir) {
  # srcDir <- "data"
  # y
  for(fpath in list.files(srcDir, pattern="_prof_y.txt$", full.names=T)) {
    cat("read ",fpath,"\n")
    # === read profile ===
    prof <- read_profiles_csv(fpath)
    g_prof <<- prof
    
    # === profile bbox ===
    flag_1line <- 1
    profiles_y_bbox(prof, sub("_prof_y.txt","_bbox_y.txt",fpath), flag_1line)
  }
  
  # x
  for(fpath in list.files(srcDir, pattern="_prof_x.txt$", full.names=T)) {
    cat("read ",fpath,"\n")
    # === read profile ===
    prof <- read_profiles_csv(fpath)
    g_prof <<- prof

    # === profile bbox ===
    flag_1line <- 1
    profiles_x_bbox(prof, sub("_prof_x.txt","_bbox_x.txt",fpath), flag_1line)
  }

  # z
  for(fpath in list.files(srcDir, pattern="_prof_z.txt$", full.names=T)) {
    cat("read ",fpath,"\n")
    # === read profile ===
    prof <- read_profiles_csv(fpath)
    g_prof <<- prof

    # === profile bbox ===
    flag_1line <- 1
    profiles_z_bbox(prof, sub("_prof_z.txt","_bbox_z.txt",fpath), flag_1line)
  }
}

slicer3d <- function(srcDir, ndiv) {
  objs_profies_xyz(srcDir, ndiv)
  profies_xyz_bbox(srcDir)
}
# slicer3d("data", 10)
