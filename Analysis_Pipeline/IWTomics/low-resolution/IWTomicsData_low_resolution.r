setClassUnion(".listOrNULL",c("list","NULL"))
setClass("IWTomicsData",
         slots=c(metadata="list",
                 regions="GRangesList",
                 alignment="character",
                 features="list",
                 length_features="list",
                 test=".listOrNULL"),
         prototype=list(metadata=list(region_datasets=data.frame(),
                                      feature_datasets=data.frame()),
                        regions=GRangesList(),
                        alignment='center',
                        features=list(),
                        length_features=list(),
                        test=NULL),
         validity=function(object){
           if(!(alignment(object) %in% c('left','right','center','scale')))
             return(paste0('invalid alignment type \'',alignment(object),'\'. Available alignments are \'left\', \'right\', \'center\' and \'scale\'.'))
           if(nRegions(object)!=0)
             if(!identical(names(object@regions),idRegions(object)))
               return('invalid regions. Names should be identical to region IDs.')
           if(nFeatures(object)!=0){
             if(!identical(names(object@features),idFeatures(object)))
               return('invalid features. Names should be identical to feature IDs.')
             if(TRUE %in% lapply(object@features,function(x) !identical(names(x),idRegions(object))))
               return('invalid features. Each feature should have names identical to region IDs.')
           }
           if(length(object@length_features)!=0){
             if(!identical(names(object@length_features),idFeatures(object)))
               return('invalid length_features. Names should be identical to feature IDs.')
             if(TRUE %in% lapply(object@length_features,function(x) !identical(names(x),idRegions(object))))
               return('invalid length_features. Each feature should have names identical to region IDs.')
           }
           if(!is.null(object@test)){
             input=testInput(object)
             if(sum(!(input$id_region1 %in% idRegions(object))))
               return('invalid test id_region1. The IDs provided are not region IDs.')
             if(!is.null(input$id_region2)){
               if(length(input$id_region2)!=length(input$id_region1))
                 return('invalid test id_region2. It must have the same length of id_region1.')
               if(sum(!(setdiff(input$id_region2,'') %in% idRegions(object))))
                 return('invalid test id_region2. The IDs provided are not region IDs.')
             }
             if(sum(!(input$id_features_subset %in% idFeatures(object))))
               return('invalid test id_features_subset. The IDs provided are not feature IDs.')
             if(!(input$statistics %in% c('mean','median','variance','quantile')))
               return('invalid \'statistics\'. Available test statistics are \'mean\', \'median\', \'variance\' and \'quantile\'.')
           }
           return(TRUE)
         })


# Constructor method
setGeneric("IWTomicsData",function(x,y,...) standardGeneric("IWTomicsData"))

setMethod("IWTomicsData",c("GRangesList","list"),
          function(x,y,alignment='center',
                   id_regions=NULL,name_regions=NULL,
                   id_features=NULL,name_features=NULL,length_features=NULL){
            if(!(alignment %in% c('left','right','center','scale')))
              stop('invalid alignment type \'',alignment,'\'. Available alignments are \'left\', \'right\', \'center\' and \'scale\'.')
            regions=x
            features=y
            if((FALSE %in% unlist(lapply(features,is.list)))|
               (FALSE %in% unlist(lapply(features,function(feature) lapply(feature,is.matrix)))))
              stop('invalid \'features\'. List of matrix lists expected.')
            if(sum(unlist(lapply(features,length))!=length(regions)))
              stop('Invalid \'features\'. A matrix for each region dataset expected.')
            if(FALSE %in% unlist(lapply(features,function(feature) length(unique(unlist(lapply(feature,nrow))))==1)))
              stop('Invalid \'features\'. Matrices for the same feature should have the same row number.')
            
            if(is.null(id_regions))
              id_regions=names(regions)
            if(is.null(id_regions))
              id_regions=paste0('rgn',seq_along(regions))
            if(is.null(name_regions))
              name_regions=id_regions
            region_datasets=data.frame(name=name_regions,file=NA,size=unlist(lapply(regions,length)),row.names=id_regions,stringsAsFactors=FALSE)
            for(id_region in id_regions){
              if(TRUE %in% duplicated(regions[[id_region]]))
                stop('duplicated regions in ',id_region,'.')
              if(!identical(disjoin(regions[[id_region]]),sort(regions[[id_region]])))
                warning('overlapping regions in ',id_region,'.')
            }
            if(is.null(id_features))
              id_features=names(features)
            if(is.null(id_features))
              id_features=paste0('ftr',seq_along(features))
            if(is.null(name_features))
              name_features=id_features
            feature_datasets=data.frame(name=name_features,row.names=id_features,stringsAsFactors=FALSE)
            feature_datasets=cbind(feature_datasets,matrix(NA,ncol=length(id_regions),nrow=length(id_features)))
            names(feature_datasets)=c('name',paste0('file_',id_regions))
            length_features_new=list()
            for(id_feature in id_features){
              length_feature=list()
              resolution=c()
              for(id_region in id_regions){
                if(is.null(length_features[[id_feature]][[id_region]])){
                  ranges=apply(features[[id_feature]][[id_region]],2,
                               function(feature){
                                 notNA=which(!is.na(feature))
                                 return(c(notNA[1],notNA[length(notNA)]))
                               })
                  length_feature[[id_region]]=ranges[2,]-ranges[1,]+1
                } else {
                  length_feature[[id_region]]=length_features[[id_feature]][[id_region]]
                }
                names(length_feature[[id_region]])=NULL
                feature_resolution=unique(width(regions[[id_region]])/length_feature[[id_region]])
                if(length(feature_resolution)>1)
                  warning('different size windows for feature \'',id_feature,'\'.')
                resolution=c(resolution,feature_resolution[1])
              }
              resolution=unique(resolution)
              if(length(resolution)>1)
                warning(paste0('Different size windows for feature \'',id_feature,'\'.'))
              feature_datasets[id_feature,'resolution']=resolution[1]
              if(max(unlist(length_feature))!=nrow(features[[id_feature]][[1]]))
                stop(paste0('Invalid features. Row of NA in the matrices of feature \'',id_feature,'\'.'))
              length_features_new[[id_feature]]=length_feature
            }
            
            new("IWTomicsData",metadata=list(region_datasets=region_datasets,feature_datasets=feature_datasets),
                regions=regions,alignment=alignment,features=features,length_features=length_features_new)
          })

setMethod("IWTomicsData",c("character","data.frame"),
          function(x,y,alignment='center',
                   id_regions=NULL,name_regions=NULL,id_features=NULL,name_features=NULL,
                   path=NULL,start.are.0based=TRUE,header=FALSE,...){
            if(!(alignment %in% c('left','right','center','scale')))
              stop('invalid alignment type \'',alignment,'\'. Available alignments are \'left\', \'right\', \'center\' and \'scale\'.')
            
            file_regions=x
            if(is.null(id_regions))
              id_regions=file_path_sans_ext(file_regions)
            if(is.null(name_regions))
              name_regions=id_regions
            region_datasets=data.frame(name=name_regions,file=file_regions,row.names=id_regions,stringsAsFactors=FALSE)
            names(name_regions)=id_regions
            
            file_features=y
            if(ncol(file_features)<length(file_regions))
              stop('invalid file_features.')
            dataset.in.files=id_regions %in% names(file_features)
            if(FALSE %in% dataset.in.files){
              if(TRUE %in% dataset.in.files)
                stop('invalid file_features.')
              names(file_features)=id_regions
            }
            if(nrow(file_features)!=1){
              file_features=data.frame(apply(file_features,2,as.character),stringsAsFactors=FALSE)
            }else{
              file_features=data.frame(t(apply(file_features,2,as.character)),stringsAsFactors=FALSE)
            }
            if(is.null(id_features))
              id_features=file_path_sans_ext(file_features[,1])
            if(is.null(name_features))
              name_features=id_features
            feature_datasets=data.frame(name=name_features,file=file_features,row.names=id_features,stringsAsFactors=FALSE)
            colnames(feature_datasets)[-1]=paste0('file_',id_regions)
            names(name_features)=id_features
            
            if(!is.null(path)){
              file_regions=file.path(path,file_regions)
              if(nrow(file_features)!=1){
                file_features=data.frame(apply(file_features,2,function(file) file.path(path,file)),stringsAsFactors=FALSE)
              }else{
                file_features=data.frame(t(apply(file_features,2,function(file) file.path(path,file))),stringsAsFactors=FALSE)
              }
            }
            names(file_regions)=id_regions
            rownames(file_features)=id_features
            
            regions=list()
            for(id_region in id_regions){
              message('Reading region dataset \'',name_regions[id_region],'\'...')
              file=file_regions[id_region]
              tmp=read.delim(file,header=header,stringsAsFactors=FALSE,...)
              if(ncol(tmp)<3)
                stop('invalid format in ',file,'.')
              if(header){
                if(sum(!(c('chr','start','end') %in% names(tmp))))
                  stop('invalid variable names in ',file,'. Variable names should be \'chr\', \'start\' and \'end\'.')
              }else{
                names(tmp)[1:3]=c('chr','start','end')
              }
              regions[[id_region]]=makeGRangesFromDataFrame(tmp,starts.in.df.are.0based=start.are.0based)
              rm(tmp)
              if(TRUE %in% duplicated(regions[[id_region]]))
                stop('duplicated regions in ',file,'.')
              if(!isDisjoint(regions[[id_region]]))
                warning('overlapping regions in ',file,'.')
            }
            region_datasets$size=unlist(lapply(regions,length))
            regions=GRangesList(regions)
            
            features=list()
            length_features=list()
            for(id_feature in id_features){
              message('Reading feature \'',name_features[id_feature],'\'...')
              feature.matrices=list()
              resolution=c()
              for(id_region in id_regions){
                message('   Region dataset \'',name_regions[id_region],'\'...')
                tmp=read.delim(file_features[id_feature,id_region],header=header,stringsAsFactors=FALSE,...)
                if(ncol(tmp)<4)
                  stop('invalid format in ',file_features[id_feature,id_region],'.')
                measure.window=ifelse(ncol(tmp)==4,TRUE,FALSE)
                if(header){
                  if(sum(!(c('chr','start','end') %in% names(tmp))))
                    stop('Invalid variable names in ',file_features[id_feature,id_region],'. Variable names should be \'chr\', \'start\' and \'end\'.')
                }else{
                  names(tmp)[1:4]=c('chr','start','end','measure')
                }
                tmp=makeGRangesFromDataFrame(tmp,starts.in.df.are.0based=start.are.0based,keep.extra.columns=TRUE)
                if(TRUE %in% duplicated(tmp))
                  stop('duplicated windows in ',file_features[id_feature,id_region],'.')
                if(!isDisjoint(tmp))
                  if(measure.window){
                    stop('overlapping windows in ',file_features[id_feature,id_region],'.')
                  }else{
                    warning('overlapping windows in ',file_features[id_feature,id_region],'.')
                  }
                if(measure.window){
                  feature.resolution=unique(width(tmp))
                }else{
                  match.region=match(regions[[id_region]],tmp)
                  if(sum(is.na(match.region)))
                    stop('not all regions in datasets ',id_region,' are present in ',file_features[id_feature,id_region],'.')
                  tmp=tmp[match.region,]
                  skip=list(...)[['skip']]
                  if(is.null(skip))
                    skip=0
                  length.tmp=(count.fields(file_features[id_feature,id_region],sep="\t",skip=header+skip)-3)[match.region]
                  feature.resolution=unique(width(regions[[id_region]])/length.tmp)
                }
                if(length(feature.resolution)>1)
                  if(measure.window){
                    warning('different size windows in ',file_features[id_feature,id_region],'.')
                  }else{
                    warning('different size windows in ',file_features[id_feature,id_region],'.')
                  }
                resolution=c(resolution,feature.resolution[1])
                if(length(setdiff(regions[[id_region]],tmp))>0)
                  warning('windows in ',file_features[id_feature,id_region],' do not cover all regions. 
                          Put NA in the file to indicate Not Available measurements.')
                if(measure.window){
                  if(region_datasets[id_region,'size']>100){
                    core.number <- min(floor(region_datasets[id_region,'size']/50),detectCores()) # at least 50 regions each node
                    n_group=max(floor(region_datasets[id_region,'size']/core.number),50)
                    groups=c(rep.int(1,region_datasets[id_region,'size']-core.number*n_group),rep(seq.int(core.number),each=n_group))
                    cl <- makeCluster(core.number)
                    feature.matrices[[id_region]]=Reduce(c,parLapply(cl,split(regions[[id_region]],groups),
                                                                     function(region,tmp) lapply(region,function(region,tmp) subsetByOverlaps(tmp,region)$measure,tmp=tmp),tmp=tmp))
                    stopCluster(cl)
                  }else{
                    feature.matrices[[id_region]]=lapply(regions[[id_region]],
                                                         function(region) subsetByOverlaps(tmp,region)$measure)
                  }
                }else{
                  feature.matrices[[id_region]]=mapply(function(tmp,length.tmp) tmp[seq_len(length.tmp)],split(as.matrix(mcols(tmp)),seq_along(length.tmp)),length.tmp,SIMPLIFY=FALSE)
                }
                names(feature.matrices[[id_region]])=NULL
                rm(tmp)
              }
              resolution=unique(resolution)
              if(length(resolution)>1)
                if(measure.window){
                  warning('Different size windows for feature \'',id_feature,'\'.')
                  resolution=min(resolution)
                }else{
                  warning('Different size windows for feature \'',id_feature,'\'.')
                  resolution=min(resolution)
                }
              feature_datasets[id_feature,'resolution']=resolution
              
              length=lapply(feature.matrices,function(feature.matrix) unlist(lapply(feature.matrix,length),use.names=FALSE))
              if(length(unique(unlist(length)))==1){
                feature.matrices=lapply(feature.matrices,function(feature) do.call(cbind,feature))
              }else{
                length.max=max(unlist(length))
                if(alignment %in% c('left','scale')){
                  length.NA.right=lapply(length,function(length) length.max-length)
                  for(id_region in id_regions)
                    feature.matrices[[id_region]]=mapply(function(feature,length.NA.right) c(feature,rep.int(NA,length.NA.right)),
                                                         feature.matrices[[id_region]],length.NA.right[[id_region]])
                }
                if(alignment=='right'){
                  length.NA.left=lapply(length,function(length) length.max-length)
                  for(id_region in id_regions)
                    feature.matrices[[id_region]]=mapply(function(feature,length.NA.left) c(rep.int(NA,length.NA.left),feature),
                                                         feature.matrices[[id_region]],length.NA.left[[id_region]])
                }
                if(alignment=='center'){
                  center=lapply(length,'%/%',2) 
                  # the alignment is approximate if there are regions with an odd number of windows and regions with an even number of regions
                  length.NA.right=lapply(mapply('-',length,center,SIMPLIFY=FALSE),function(lenght_center) (length.max-length.max%/%2)-lenght_center)
                  length.NA.left=lapply(center,function(center) length.max%/%2-center)
                  for(id_region in id_regions)
                    feature.matrices[[id_region]]=mapply(function(feature,length.NA.left,length.NA.right) as.matrix(c(rep.int(NA,length.NA.left),feature,rep.int(NA,length.NA.right))),
                                                         feature.matrices[[id_region]],length.NA.left[[id_region]],length.NA.right[[id_region]])
                }
              }
              
              features[[id_feature]]=feature.matrices
              length_features[[id_feature]]=length
          }
            
            new("IWTomicsData",metadata=list(region_datasets=region_datasets,feature_datasets=feature_datasets),
                regions=regions,alignment=alignment,features=features,length_features=length_features)
            })

setMethod("IWTomicsData",c("character","character"),
          function(x,y,...){
            y=data.frame(matrix(y,ncol=length(x),nrow=length(y)),stringsAsFactors=FALSE)
            IWTomicsData(x,y,...)
          })

setMethod("IWTomicsData",c("character","matrix"),
          function(x,y,...){
            y=data.frame(y,stringsAsFactors=FALSE)
            IWTomicsData(x,y,...)
          })

# Accessors
setGeneric("nRegions",function(x,...) standardGeneric("nRegions"))
setMethod("nRegions","IWTomicsData",function(x) length(x@regions))

setGeneric("nFeatures",function(x,...) standardGeneric("nFeatures"))
setMethod("nFeatures","IWTomicsData",function(x) length(x@features))

setMethod("dim","IWTomicsData",function(x) c(nRegions(x),nFeatures(x)))

setGeneric("lengthRegions",function(x,...) standardGeneric("lengthRegions"))
setMethod("lengthRegions","IWTomicsData",function(x) unlist(lapply(x@regions,length)))

setGeneric("lengthFeatures",function(x,...) standardGeneric("lengthFeatures"))
setMethod("lengthFeatures","IWTomicsData",function(x) x@length_features)

setGeneric("resolution",function(x,...) standardGeneric("resolution"))
setMethod("resolution","IWTomicsData",
          function(x){
            res=metadata(x)$feature_datasets$resolution
            names(res)=idFeatures(x)
            return(res)
          })

setMethod("metadata","IWTomicsData",function(x) x@metadata)

setGeneric("regions",function(x,...) standardGeneric("regions"))
setMethod("regions","IWTomicsData",function(x) x@regions)

setGeneric("features",function(x,...) standardGeneric("features"))
setMethod("features","IWTomicsData",function(x) x@features)

setGeneric("idRegions",function(x,...) standardGeneric("idRegions"))
setMethod("idRegions","IWTomicsData",function(x) row.names(metadata(x)$region_datasets))

setGeneric("idFeatures",function(x,...) standardGeneric("idFeatures"))
setMethod("idFeatures","IWTomicsData",function(x) row.names(metadata(x)$feature_datasets))

setGeneric("nameRegions",function(x,...) standardGeneric("nameRegions"))
setMethod("nameRegions","IWTomicsData",
          function(x){
            name=metadata(x)$region_datasets$name
            names(name)=idRegions(x)
            name
          })

setGeneric("nameFeatures",function(x,...) standardGeneric("nameFeatures"))
setMethod("nameFeatures","IWTomicsData",
          function(x){
            name=metadata(x)$feature_datasets$name
            names(name)=idFeatures(x)
            name
          })

setGeneric("alignment",function(x,...) standardGeneric("alignment"))
setMethod("alignment","IWTomicsData",function(x) x@alignment)

setGeneric("testInput",function(x,...) standardGeneric("testInput"))
setMethod("testInput","IWTomicsData",function(x) x@test$input)

setGeneric("nTests",function(x,...) standardGeneric("nTests"))
setMethod("nTests","IWTomicsData",function(x) length(testInput(x)$id_region1))

setGeneric("idRegionsTest",function(x,test,...) standardGeneric("idRegionsTest"))
setMethod("idRegionsTest",c("IWTomicsData","vector"),
          function(x,test){
            if(nTests(x)==0)
              return(NULL)
            if(sum(!(test %in% 1:nTests(x))))
              stop('invalid test number.')
            id=lapply(test,
                      function(i) c(testInput(x)$id_region1[i],
                                    ifelse(is.null(testInput(x)$id_region2)||(testInput(x)$id_region2[i]==''),'',testInput(x)$id_region2[i])))
            names(id)=paste0('test',test)
            return(id)
          })
setMethod("idRegionsTest",c("IWTomicsData","missing"),
          function(x){
            if(nTests(x)==0)
              return(NULL)
            idRegionsTest(x,1:nTests(x))
          })

setGeneric("idFeaturesTest",function(x,...) standardGeneric("idFeaturesTest"))
setMethod("idFeaturesTest","IWTomicsData",function(x) testInput(x)$id_features_subset)

setGeneric("adjusted_pval",function(x,test,id_features_subset,scale_threshold,...) standardGeneric("adjusted_pval"))
setMethod("adjusted_pval",c("IWTomicsData","vector","character","vector"),
          function(x,test,id_features_subset,scale_threshold){
            if(nTests(x)==0)
              return(NULL)
            if(sum(!(test %in% 1:nTests(x))))
              stop('invalid test number.')
            if(sum(!(id_features_subset %in% idFeaturesTest(x))))
              stop('invalid id_features_subset.')
            if(is.list(scale_threshold)){
              scale_threshold=lapply(scale_threshold,
                                     function(scale){
                                       scale=as.list(rep(scale,length.out=length(id_features_subset)))
                                       names(scale)=id_features_subset
                                       return(scale)
                                     })
            }else{
              scale_threshold=lapply(test,
                                     function(i){
                                       scale=as.list(rep(scale_threshold,length.out=length(id_features_subset)))
                                       names(scale)=id_features_subset
                                       return(scale)
                                     })
            }
            pval=mapply(function(results,scale) mapply(function(result,scale){
              pval=result$adjusted_pval_matrix
              if((scale<1)||(result$max_scale<scale)){
                warning('invalid scale_threshold. Setting it to the default value.',call.=FALSE,immediate.=TRUE)
                scale=result$max_scale
              }
              pval=pval[ncol(pval)-scale+1,]
              return(pval)
            },results,scale,SIMPLIFY=FALSE),
            .testResults(x,test,id_features_subset),scale_threshold,SIMPLIFY=FALSE)
            names(pval)=paste0('test',test)
            return(pval)
          })
setMethod("adjusted_pval",c("IWTomicsData","missing","character","vector"),
          function(x,test,id_features_subset,scale_threshold){
            if(nTests(x)==0)
              return(NULL)
            adjusted_pval(x,1:nTests(x),id_features_subset,scale_threshold)
          })
setMethod("adjusted_pval",c("IWTomicsData","vector","missing","vector"),
          function(x,test,id_features_subset,scale_threshold){
            if(nTests(x)==0)
              return(NULL)
            adjusted_pval(x,test,idFeaturesTest(x),scale_threshold)
          })
setMethod("adjusted_pval",c("IWTomicsData","missing","missing","vector"),
          function(x,test,id_features_subset,scale_threshold){
            if(nTests(x)==0)
              return(NULL)
            adjusted_pval(x,1:nTests(x),idFeaturesTest(x),scale_threshold)
          })
setMethod("adjusted_pval",c("IWTomicsData","vector","character","missing"),
          function(x,test,id_features_subset){
            if(nTests(x)==0)
              return(NULL)
            if(sum(!(test %in% 1:nTests(x))))
              stop('invalid test number.')
            if(sum(!(id_features_subset %in% idFeaturesTest(x))))
              stop('invalid id_features_subset.')
            pval=lapply(x@test$result[test],function(results) lapply(results[id_features_subset],function(result) result$adjusted_pval))
            names(pval)=paste0('test',test)
            return(pval)
          })
setMethod("adjusted_pval",c("IWTomicsData","missing","character","missing"),
          function(x,test,id_features_subset){
            if(nTests(x)==0)
              return(NULL)
            adjusted_pval(x,1:nTests(x),id_features_subset)
          })
setMethod("adjusted_pval",c("IWTomicsData","vector","missing","missing"),
          function(x,test){
            if(nTests(x)==0)
              return(NULL)
            adjusted_pval(x,test,idFeaturesTest(x))
          })
setMethod("adjusted_pval",c("IWTomicsData","missing","missing","missing"),
          function(x){
            if(nTests(x)==0)
              return(NULL)
            adjusted_pval(x,1:nTests(x),idFeaturesTest(x))
          })

setGeneric(".testResults",function(x,test,id_features_subset,...) standardGeneric(".testResults"))
setMethod(".testResults",c("IWTomicsData","vector","character"),
          function(x,test,id_features_subset){
            if(nTests(x)==0)
              return(NULL)
            if(sum(!(test %in% 1:nTests(x))))
              stop('invalid test number.')
            if(sum(!(id_features_subset %in% idFeaturesTest(x))))
              stop('invalid id_features_subset.')
            lapply(x@test$result[test],function(results) results[id_features_subset])
          })
setMethod(".testResults",c("IWTomicsData","missing","character"),
          function(x,test,id_features_subset){
            if(nTests(x)==0)
              return(NULL)
            .testResults(x,1:nTests(x),id_features_subset)
          })
setMethod(".testResults",c("IWTomicsData","vector","missing"),
          function(x,test,id_features_subset){
            if(nTests(x)==0)
              return(NULL)
            .testResults(x,test,idFeaturesTest(x))
          })
setMethod(".testResults",c("IWTomicsData","missing","missing"),
          function(x){
            if(nTests(x)==0)
              return(NULL)
            .testResults(x,1:nTests(x),idFeaturesTest(x))
          })


# Subset method
setMethod("[",c("IWTomicsData","ANY","ANY","ANY"),
          function(x,i,j,...,drop=TRUE){
            if(!missing(i)){
              if(is.numeric(i)&(FALSE %in% (i %in% seq_len(nRegions(x)))))
                stop('undefined regions selected')
              if(is.character(i)&(FALSE %in% (i %in% idRegions(x))))
                stop('undefined regions selected')
            }
            if(!missing(j)){
              if(is.numeric(j)&(FALSE %in% (j %in% seq_len(nFeatures(x)))))
                stop('undefined features selected')
              if(is.character(j)&(FALSE %in% (j %in% idFeatures(x))))
                stop('undefined features selected')
            }
            region_datasets=x@metadata$region_datasets[i,]
            feature_datasets=x@metadata$feature_datasets[j,c('name',paste0('file_',row.names(region_datasets)),'resolution')]
            regions_new=x@regions[row.names(region_datasets)]
            length_features=lapply(x@length_features[row.names(feature_datasets)],function(length) length[row.names(region_datasets)])
            length_features_max=lapply(length_features,function(length) max(unlist(length)))
            features_new=mapply(function(feature,length){
              feature=feature[row.names(region_datasets)]
              if(alignment(x) %in% c('left','scale'))
                return(lapply(feature,function(M) as.matrix(M[seq_len(length),])))
              if(alignment(x)=='right')
                return(lapply(feature,function(M) as.matrix(M[seq_len(length)+nrow(M)-length,])))
              if(alignment(x)=='center')
                return(lapply(feature,function(M) as.matrix(M[seq_len(length)+nrow(M)%/%2-length%/%2,])))
            },x@features[row.names(feature_datasets)],length_features_max,SIMPLIFY=FALSE)
            initialize(x,metadata=list(region_datasets=region_datasets,feature_datasets=feature_datasets),
                       regions=regions_new,features=features_new,length_features=length_features,test=NULL)
          })


# Combine methods
setMethod("c","IWTomicsData",
          function(x,...){
            elements=list(x,...)
            if(length(elements)>2){
              c(elements[[1]],c(...))
            }else{
              x=elements[[1]]
              y=elements[[2]]
              alignment_new=unique(alignment(x),alignment(y))
              if(length(alignment_new)>1)
                stop('merging not possible, different types of alignment present.')
              regions_id=union(idRegions(x),idRegions(y))
              features_id=union(idFeatures(x),idFeatures(y))
              
              regions_common=intersect(idRegions(x),idRegions(y))
              for(region_id in regions_common){
                if(length(setdiff(regions(x)[[region_id]],regions(y)[[region_id]])))
                  stop('merging not possible, region dataset \'',region_id,'\' differs in the IWTomicsData objects.')
                overlaps=as.matrix(findOverlaps(regions(x)[[region_id]],regions(y)[[region_id]]))
                if(sum(overlaps[,1]!=1:length(regions(x)[[region_id]])))
                  stop('merging not possible, region dataset \'',region_id,'\' differs in the IWTomicsData objects.')
                if(sum(overlaps[,2]!=1:length(regions(y)[[region_id]]))){
                  y@regions[[region_id]]=y@regions[[region_id]][overlaps[,2]]
                  for(feature_id in idFeatures(y)){
                    y@features[[feature_id]][[region_id]]=y@features[[feature_id]][[region_id]][,overlaps[,2]]
                    y@length_features[[feature_id]][[region_id]]=y@length_features[[feature_id]][[region_id]][overlaps[,2]]
                  }
                }
              }
              features_common=intersect(idFeatures(x),idFeatures(y))
              for(feature_id in features_common){
                if(metadata(x)$feature_datasets[feature_id,'resolution']!=metadata(y)$feature_datasets[feature_id,'resolution'])
                  stop('merging not possible, feature \'',feature_id,'\' resolution differs in the IWTomicsData objects.')
                for(region_id in regions_common){
                  if(!identical(features(x[region_id,feature_id]),features(y[region_id,feature_id])))
                    stop('merging not possible, feature \'',feature_id,'\' in region dataset \'',region_id,'\' differs in the IWTomicsData objects.')
                }
              }
              for(region_id in regions_id){
                features_id_present=c()
                if(region_id %in% idRegions(x))
                  features_id_present=idFeatures(x)
                if(region_id %in% idRegions(y))
                  features_id_present=c(features_id_present,idFeatures(y))
                if(!isEmpty(setdiff(features_id,features_id_present)))
                  stop('merging not possible, not all features are present for region dataset \'',region_id,'\'.')
              }
              if(isEmpty(setdiff(idRegions(x),regions_common))&isEmpty(setdiff(idFeatures(x),features_common)))
                return(y)
              if(isEmpty(setdiff(idRegions(y),regions_common))&isEmpty(setdiff(idFeatures(y),features_common)))
                return(x)
              
              region_datasets=rbind(x@metadata$region_datasets,y@metadata$region_datasets[setdiff(idRegions(y),regions_common),])
              feature_datasets=as.data.frame(matrix(NA,nrow=length(features_id),ncol=length(regions_id)+2))
              row.names(feature_datasets)=features_id
              colnames(feature_datasets)=c('name',paste0('file_',regions_id),'resolution')
              feature_datasets[idFeatures(x),c('name',paste0('file_',idRegions(x)),'resolution')]=x@metadata$feature_datasets[idFeatures(x),c('name',paste0('file_',idRegions(x)),'resolution')]
              feature_datasets[idFeatures(y),c('name',paste0('file_',idRegions(y)),'resolution')]=y@metadata$feature_datasets[idFeatures(y),c('name',paste0('file_',idRegions(y)),'resolution')]
              regions_new=c(regions(x),regions(y)[setdiff(idRegions(y),regions_common)])
              features_new=lapply(features_id,
                                  function(feature_id){
                                    feature=vector('list',length(regions_id))
                                    names(feature)=regions_id
                                    return(feature)
                                  })
              names(features_new)=features_id
              length_features=features_new
              for(feature_id in idFeatures(x))
                for(region_id in idRegions(x)){
                  features_new[[feature_id]][[region_id]]=x@features[[feature_id]][[region_id]]
                  length_features[[feature_id]][[region_id]]=x@length_features[[feature_id]][[region_id]]
                }
              for(feature_id in idFeatures(y))
                for(region_id in idRegions(y)){
                  features_new[[feature_id]][[region_id]]=y@features[[feature_id]][[region_id]]
                  length_features[[feature_id]][[region_id]]=y@length_features[[feature_id]][[region_id]]
                }
              for(feature_id in features_id){
                length=lapply(features_new[[feature_id]],nrow)
                if(length(unique(unlist(length)))!=1){
                  length.max=max(unlist(length))
                  if(alignment_new %in% c('left','scale')){
                    length.NA.right=lapply(length,function(length) length.max-length)
                    for(region_id in regions_id)
                      if(length(length.NA.right[[region_id]])>0)
                        features_new[[feature_id]][[region_id]]=rbind(features_new[[feature_id]][[region_id]],
                                                                      matrix(NA,nrow=length.NA.right[[region_id]],ncol=ncol(features_new[[feature_id]][[region_id]])))
                  }
                  if(alignment_new=='right'){
                    length.NA.left=lapply(length,function(length) length.max-length)
                    for(region_id in regions_id)
                      if(length(length.NA.left[[region_id]])>0)
                        features_new[[feature_id]][[region_id]]=rbind(matrix(NA,nrow=length.NA.left[[region_id]],ncol=ncol(features_new[[feature_id]][[region_id]])),
                                                                      features_new[[feature_id]][[region_id]])
                  }
                  if(alignment_new=='center'){
                    center=lapply(length,'%/%',2) 
                    # the alignment is approximate if there are regions with an odd number of windows and regions with an even number of regions
                    length.NA.right=lapply(mapply('-',length,center,SIMPLIFY=FALSE),function(lenght_center) (length.max-length.max%/%2)-lenght_center)
                    length.NA.left=lapply(center,function(center) length.max%/%2-center)
                    for(region_id in regions_id)
                      if(length(length.NA.right[[region_id]])>0)
                        features_new[[feature_id]][[region_id]]=rbind(matrix(NA,nrow=length.NA.left[[region_id]],ncol=ncol(features_new[[feature_id]][[region_id]])),
                                                                      features_new[[feature_id]][[region_id]],
                                                                      matrix(NA,nrow=length.NA.right[[region_id]],ncol=ncol(features_new[[feature_id]][[region_id]])))
                  }
                }
              }
              new("IWTomicsData",metadata=list(region_datasets=region_datasets,feature_datasets=feature_datasets),
                  regions=regions_new,alignment=alignment_new,features=features_new,length_features=length_features,test=NULL)
            }
          })

setMethod("merge",c("IWTomicsData","IWTomicsData"),function(x,y,...) c(x,y,...))

setMethod("rbind","IWTomicsData",
          function(...){
            elements=list(...)
            features_id=lapply(elements,idFeatures)
            equal.features_id=unlist(lapply(seq_along(features_id)[-1],function(i) identical(features_id[[i-1]],features_id[[i]])))
            if(FALSE %in% equal.features_id)
              stop('merging not possible, features differs in the IWTomicsData objects.')
            c(...)})

setMethod("cbind","IWTomicsData",
          function(...){
            elements=list(...)
            regions_id=lapply(elements,idRegions)
            equal.regions_id=unlist(lapply(seq_along(regions_id)[-1],function(i) identical(regions_id[[i-1]],regions_id[[i]])))
            if(FALSE %in% equal.regions_id)
              stop('merging not possible, locations differs in the IWTomicsData objects.')
            c(...)})
