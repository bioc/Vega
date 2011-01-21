plotSegmentation <- function(CNVdata, segmentation, chromosomes, opt=0){

# opt=0 plots the means
# opt=1 plots the CN

	num_of_chroms <- length(chromosomes);
	
	if(num_of_chroms<=6){
		rows <- 1;
		cols <- num_of_chroms;
	}else{
		cols <- 6;
		tmp1 <- round(num_of_chroms/cols);
		tmp2 <- num_of_chroms/cols;
		diff_tmp <- tmp1-tmp2;
		
		if(diff_tmp==0){
			rows <- tmp1;
		}else{
			rows <- tmp1+1;
		}
	}
	
	par(mfrow=c(rows,cols));

	for(i in 1:num_of_chroms){
		curr_chr <- chromosomes[i];
		
		data_table <- CNVdata[which(CNVdata[,1]==curr_chr),];
		if(class(data_table)!="matrix"){
			data_table <- t(as.matrix(data_table));
		}

		curr_seg <- segmentation[which(segmentation[,1]==curr_chr),];		
		if(class(curr_seg)!="matrix"){
			curr_seg <- t(as.matrix(curr_seg));
		}
		
		min_lrr <- min(as.numeric(data_table[,4]))-0.5;
		max_lrr <- max(as.numeric(data_table[,4]))+0.5;
		if(opt==1 && min_lrr > -1.5){
			min_lrr <- -1.5;
		}
		if(opt==1 && max_lrr < 1.5){
			max_lrr <- 1.5;
		}
	plot(as.numeric(data_table[,4]), col="blue", xaxt="n", xlab="", ylab="LRR", main=paste("chromosome", curr_chr), ylim=c(min_lrr, max_lrr));
		chr_size <- sum(as.numeric(curr_seg[,4]));
		tmp <- 0*c(1:chr_size);
		curr = 1;
		for(i in 1:nrow(curr_seg)){
			if(opt==0){
				label = as.numeric(curr_seg[i,5]);
			}
			if(opt==1){
				label = as.numeric(curr_seg[i,6]);
			}
			tmp[curr:(curr+as.numeric(curr_seg[i,4]))] <- label;
			curr = (curr+as.numeric(curr_seg[i,4]));
		}
		
		
		lines(c(1:length(tmp)), tmp, col="red", lwd=2);
	}
	
	
}
