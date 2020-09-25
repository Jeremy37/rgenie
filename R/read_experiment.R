
get_region_reads_from_bam = function(name, replicate, bam_file, chr, start, end, ref_sequence, min_mapq) {
  region = IRanges::IRangesList(seq1 = IRanges::IRanges(start, end))
  names(region) = chr
  read_params = Rsamtools::ScanBamParam(which = region,
                                         what = c("rname", "strand", "pos", "qwidth", "cigar", "seq"),
                                         mapqFilter = min_mapq,
                                         flag = Rsamtools::scanBamFlag(isSecondaryAlignment = FALSE, isUnmappedQuery = FALSE))
  bam_data = Rsamtools::scanBam(bam_file, param = read_params)
  return(bam_data[[1]])
}

# this function takes a .sam read and registers to the reference genome, relative
# to the region of interest in the reference
# '-' not sequenced nucleotides
# '*' deleted nucleotides
# if the cigars contains H or I the variable 'read_ok' will be False
# the returned read should be the same length as the reference genome
register_read = function(relative_pos, sam_cigar, sam_read, region_length) {
  sam_cigar_parsed = parse_cigar(sam_cigar)
  read_ok = TRUE
  read_status = ""
  expanded_cigar = ""
  cigar_bits = ""
  cursor_read = 1
  match_length = 0
  del_length = 0
  read_span_to = relative_pos + 1
  left_softclip_len = 0
  right_softclip_pos = -1
  right_softclip_len = 0

  for (i in 1:length(sam_cigar_parsed$letters)) {
    type = sam_cigar_parsed$letters[i]

    if (type == 'M') {
      match_length = match_length + sam_cigar_parsed$numbers[i]
      cigar_bits[i] = substr(sam_read, cursor_read, cursor_read + sam_cigar_parsed$numbers[i] - 1)
      #expanded_cigar = paste0(expanded_cigar, substr(sam_read, cursor_read, cursor_read + sam_cigar_parsed$numbers[i] - 1))
      cursor_read = cursor_read + sam_cigar_parsed$numbers[i]
      read_span_to = read_span_to + sam_cigar_parsed$numbers[i]
    }
    else if (type == 'D') {
      cigar_bits[i] = strrep('*', sam_cigar_parsed$numbers[i])
      #expanded_cigar = paste0(expanded_cigar, strrep('*', sam_cigar_parsed$numbers[i]))
      del_length = del_length + sam_cigar_parsed$numbers[i]
      read_span_to = read_span_to + sam_cigar_parsed$numbers[i]
    } else if (type == 'N') {
      read_status = "spliced"
      cigar_bits[i] = strrep('-', sam_cigar_parsed$numbers[i])
    } else if (type == 'S') {
      if (left_softclip_len > 0) {
        right_softclip_pos = read_span_to
        right_softclip_len = sam_cigar_parsed$numbers[i]
      } else {
        left_softclip_len = sam_cigar_parsed$numbers[i]
      }
      cursor_read = cursor_read + sam_cigar_parsed$numbers[i]
      read_status = "softclipped"
    } else if (type == 'H') {
      read_status = "hardclipped"
      read_ok = FALSE
      break
    } else if (type == 'I') {
      read_status = "insertion"
      read_ok = FALSE
      break
    } else {
      read_ok = FALSE
      break
    }
  }
  expanded_cigar = paste0(cigar_bits, collapse="")
  registered_read = ""
  if (read_ok) {
    if (relative_pos < 0) {
      if ((nchar(expanded_cigar) - relative_pos) > 0) {
        registered_read = substr(expanded_cigar, (-relative_pos)+1, nchar(expanded_cigar))
      }
    } else if (relative_pos == 0) {
      registered_read = expanded_cigar
    } else {
      registered_read = paste0(strrep('-', relative_pos), expanded_cigar)
    }

    len = nchar(registered_read)
    if (len > region_length) {
      registered_read = substr(registered_read, 1, region_length)
    } else if (len < region_length){
      registered_read = paste0(registered_read, strrep('-', region_length - len))
    }
  }

  registered_left_softclip = vector(mode = "integer", length = region_length)
  registered_right_softclip = vector(mode = "integer", length = region_length)
  if (left_softclip_len > 0 & relative_pos > 0) {
    # relative_pos is the read start position relative to the ref sequence, but soft-clipping
    # occurs to the left of the read start. E.g. relative_pos == 0 means the non-clipped read
    # sequence starts at the first base of the ref sequence.
    startpos = relative_pos + 1 - left_softclip_len
    if (startpos < 1) {
      left_softclip_len = left_softclip_len - (1 - startpos)
      startpos = 1
    }
    if (startpos + left_softclip_len - 1 > region_length) {
      left_softclip_len = region_length - startpos + 1
    }
    registered_left_softclip[startpos:(startpos + left_softclip_len - 1)] = 1
  }
  if (right_softclip_len > 0 & relative_pos < region_length) {
    startpos = right_softclip_pos
    if (startpos < 1) {
      right_softclip_len = right_softclip_len - (1 - startpos)
      startpos = 1
    }
    if (startpos + right_softclip_len - 1 > region_length) {
      right_softclip_len = region_length - startpos + 1
    }
    registered_right_softclip[startpos:(startpos + right_softclip_len - 1)] = 1
  }

  return(list("read" = registered_read, "left_softclip" = registered_left_softclip, "right_softclip" = registered_right_softclip,
              "read_ok" = read_ok, "read_status" = read_status, "match_length" = match_length, "deletion_length" = del_length))
}

# this function takes a string cigar "10M20D7M" and converts it to arrays of letters ['M','D','M'] and numbers [10,20,7]
parse_cigar = function(input_cigar) {
  #print(input_cigar)
  X = strsplit(input_cigar, split="")[[1]] # splits the string into array of char
  numbers = numeric()
  letters = character()
  idx_last_letter = 0
  idx_current = 1

  for (x in X) {
    if (is.na(strtoi(x))) {
      letters = c(letters, x)
      idx_1 = idx_last_letter + 1
      idx_2 = idx_current - 1
      #numbers = c(numbers, as.integer(paste0(X[idx_1:idx_2], collapse="")))
      numbers = c(numbers, strtoi(substr(input_cigar, idx_1, idx_2)))
      idx_last_letter = idx_current
    }
    idx_current = idx_current + 1
  }
  return(list("letters" = letters, "numbers" = numbers))
}


# Goes through each SAM file read and gets its sequence with respect to the reference
# zone of interest, where "-" indicates that the read does not cover the position, and
# * indicates a deletion.
get_aligned_reads = function(read_data, ref_sequence, start, end) {
  # read_data: list of results from Rsamtools::scanBam, including cigar strings and read seqs
  # ref_sequence: sequence of the amplicon
  # start, end: coordinates of the region of interest in the reference
  num_reads = length(read_data$seq)
  ref_seq_length = nchar(ref_sequence) # size of the ref genome
  if (ref_seq_length != (end - start + 1)) {
    stop(sprintf("Error: reference sequence length (%d) should be the same as the region size (%d = end - start + 1)",
                 ref_seq_length, end - start + 1))
  }

  discard = vector(mode = "logical", length = num_reads) # Preallocate a vector of discarded reads
  num_outsidewindow = 0
  for (i in 1:num_reads) {
    if (read_data$cigar[i] == "*") {
      discard[i] = TRUE
      next # Ignore this read, as it doesn't have a valid CIGAR
    }
    relative_pos = read_data$pos[i] - start
    if (abs(relative_pos) > 1000) {
      num_outsidewindow = num_outsidewindow + 1
      stop("Likely ERROR: read position and reference coordinates differ by more than 1000 bp")
    }
  }
  reads.df = tibble(seq = as.character(read_data$seq),
                    sam_position = read_data$pos,
                    sam_cigar = read_data$cigar,
                    discard = discard)
  reads.summary.df = reads.df %>% group_by(sam_position, sam_cigar, seq) %>%
    summarise(count = n(), discard = first(discard)) %>%
    dplyr::mutate(relative_pos = sam_position - start)

  discardedReads.df = reads.summary.df %>% dplyr::filter(discard == TRUE)
  keptReads.df = reads.summary.df %>% dplyr::filter(!discard)

  alignedReads = vector(mode = "character", length = nrow(keptReads.df)) # Preallocate a vector of aligned reads

  num_softclipped = 0
  num_hardclipped = 0
  num_insertion = 0

  # As we go through the aligned reads we check for possible primer dimers
  num_primerdimer = 0
  primer_dimer_start = substr(ref_sequence, 1, 20)
  primer_dimer_end = substr(ref_sequence, nchar(ref_sequence) - 19, nchar(ref_sequence))
  isPrimerDimer = function(reg_read, match_length) {
    (match_length < 40 & (substr(reg_read, 1, 20) == primer_dimer_start | substr(reg_read, nchar(reg_read) - 19, nchar(reg_read)) == primer_dimer_end))
  }

  for (i in 1:nrow(keptReads.df)) {
    # computes the read sequence within the coords of interest
    output = register_read(keptReads.df$relative_pos[i], keptReads.df$sam_cigar[i], keptReads.df$seq[i], ref_seq_length)
    if (output$read_ok) {
      alignedReads[i] = output$read
    }
    if (output$read_status == "hardclipped") {
      num_hardclipped = num_hardclipped + keptReads.df$count[i]
      if (isPrimerDimer(output$read, output$match_length)) {
        num_primerdimer = num_primerdimer + keptReads.df$count[i]
      }
    } else if (output$read_status == "insertion") {
      num_insertion = num_insertion + keptReads.df$count[i]
    } else if (output$read_status == "softclipped") {
      num_softclipped = num_softclipped + keptReads.df$count[i]
      if (isPrimerDimer(output$read, output$match_length)) {
        num_primerdimer = num_primerdimer + keptReads.df$count[i]
      }
    }
  }
  keptReads.df$region_read = alignedReads
  keptReads.df = keptReads.df %>% dplyr::filter(region_read != "") %>%
    dplyr::select(region_read, count, sam_position, sam_cigar, seq)
  # Some reads may have a distinct sequence yet have the same registered
  # read sequence, because the soft-clipped sequence could differ.

  return(list("alignedReads" = keptReads.df,
              "discardedReads" = discardedReads.df,
              "num_reads" = num_reads,
              "num_hardclipped" = num_hardclipped,
              "num_insertion" = num_insertion,
              "num_softclipped" = num_softclipped,
              "num_primerdimer" = num_primerdimer,
              "num_outsidewindow" = num_outsidewindow))
}

