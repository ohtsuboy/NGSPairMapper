//
//  main.m
// NGSPairMapTool
//
//  Created by Yoshiyuki Otsubo on 2017/09/28.
//  © 2017 Yoshiyuki Otsubo. All rights reserved.
//

// Usage Instructions
//
//    This code processes paired-end read data generated using the TruSeq Standard mRNA Sample Prep kit.
//    The reads are mapped to a reference genome, and the results are analyzed to assign read pairs to their corresponding genomic locations.
//    Unlike typical mapping tools that focus on exact mapping positions, this tool counts the entire transcript encompassing the read pair, capturing the full region spanned by the reads.
//    The primary steps include:
//
//    Prepare Reference Genome: Convert the reference genome sequence into a dictionary of k-mers.
//    Extend Sequences: Extend sequences for circular genomes.
//    Process Paired-End Reads: Read the paired-end reads, map them to the reference genome, and validate the pairs.
//    Analyze Mapped Reads: Calculate coverage depth and output results.
//
//How to Use
//
//    Set the Paths: Modify the file paths for the reference genome, read pairs, and output directory.
//    Compile the Code: Use Xcode or the command line to compile the code.
//    Run the Program: Execute the compiled program. The results will be stored in the specified output directory.
//
//Mapping Algorithm
//
//    Exact Matching of First k-mer: The algorithm maps read pairs by exact matching of the first k-mer (21-mer by default) from both reads.
//    Coverage Calculation: For each valid pair, the coverage is calculated from the 5' end of read 1 to the 5' end of read 2.
//    Distance Filtering: Only pairs with an end-to-end distance (EED) between 20 and 1200 bases are considered valid.
//    Shortest EED Selection: If multiple mapping locations are found, the pair with the shortest EED is selected.
//
//Additional Information
//
//    For each replicon sequence in the genome, the program outputs position and depth data for both strands.
//    This data can be used with other software, such as ShortReadManager (https://apps.apple.com/jp/app/shortreadmanager/id1568857700),
//    to create gene-specific expression tables or visualized using the GeneDraw tool included with the GenomeMatcher software (https://www.ige.tohoku.ac.jp/joho/portalsite/files/GenomeMatcher3.php).

#import <Foundation/Foundation.h>
#import "PairData.h"

const unichar t_table[129] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 0-9
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 10-19
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 20-29
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 30-39
    0, 0, '*', 0, 0, '-', 0, 0, 0, 0, // 40-49
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 50-59
    0, 0, 0, '?', 0, 'T', 'V', 'G', 'H', 0, // 60-69
    0, 'C', 'D', 0, 0, 'M', 0, 'K', 'N', 0, // 70-79
    0, 0, 'Y', 'S', 'A', 'A', 'B', 'W', 0, 'R', // 80-89
    0, 0, 0, 0, 0, 0, 0, 't', 'v', 'g', // 90-99
    'h', 0, 0, 'c', 'd', 0, 0, 'm', 0, 'k', // 100-109
    'n', 0, 0, 0, 'y', 's', 'a', 'a', 'b', 'w', // 110-119
    0, 'r', 0, 0, 0, 0, 0, 0, 0 // 120-128
};

NSMutableDictionary *dictionaries;
NSMutableDictionary *headerAndSeq;

NSString *complement(NSString *DNA);
void dictionayOfHeadderAndSeq(NSString *string, NSMutableDictionary *dictionary, int mode);
NSMutableDictionary* extendSequence(NSMutableDictionary *headerAndSeq, int fixedExtension);
PairData *validatePair2(NSString *st1, NSString *st2);
int mer;

NSMutableDictionary* extendSequenceIfNeeded(NSMutableDictionary *headerAndSeq) {
    NSMutableDictionary* seqs_extendedIfNeeded = [NSMutableDictionary dictionary];
    for (NSString *seqName in headerAndSeq) {
        @autoreleasepool {
            BOOL isCircular = NO;
            if ([seqName rangeOfString:@"topology=circular"].location != NSNotFound) {
                // Treat as circular if header contains this string
                isCircular = YES;
            }
            NSString *seq;
            if (isCircular) {
                NSString *s = [headerAndSeq objectForKey:seqName];
                NSString *app = [s substringToIndex:mer-1];
                seq = [NSString stringWithFormat:@"%@%@",s,app];
            } else {
                seq = [headerAndSeq objectForKey:seqName];
            }
            [seqs_extendedIfNeeded setObject:seq forKey:seqName];
        }
    }
    return seqs_extendedIfNeeded;
}

NSMutableDictionary* extendSequence(NSMutableDictionary *headerAndSeq, int fixedExtension) {
    NSMutableDictionary* superExtendSequence = [NSMutableDictionary dictionary];
    for (NSString *seqName in headerAndSeq) {
        @autoreleasepool {
            BOOL isCircular = NO;
            if ([seqName rangeOfString:@"topology=circular"].location != NSNotFound) {
                // Treat as circular if header contains this string
                isCircular = YES;
            }
            NSString *seq;
            if (isCircular) {
                NSString *s = [headerAndSeq objectForKey:seqName];
                NSString *app = [s substringToIndex:fixedExtension];
                seq = [NSString stringWithFormat:@"%@%@",s,app];
            } else {
                seq = [headerAndSeq objectForKey:seqName];
            }
            [superExtendSequence setObject:seq forKey:seqName];
        }
    }
    return superExtendSequence;
}

NSMutableDictionary* getReadyForDictionary(NSMutableDictionary *headerAndSeq) {
    dictionaries = [NSMutableDictionary dictionary];
    for (NSString *seqName in headerAndSeq) {
        @autoreleasepool {
            NSString *seq = [headerAndSeq objectForKey:seqName];
            NSMutableDictionary *dicF = [NSMutableDictionary dictionary];
            NSMutableDictionary *dicR = [NSMutableDictionary dictionary];
            [dictionaries setObject:[NSArray arrayWithObjects:dicF,dicR,nil] forKey:seqName];
            // Fragment sequence into kmers and store in dictionaries
            for (int xxx = 0 ; xxx < [seq length] - mer ; xxx++) {
                NSString *kmer = [seq substringWithRange:NSMakeRange(xxx, mer)];
                id array;
                if ((array = [dicF objectForKey:kmer])) {
                    // If kmer already exists, add position
                    [array addObject:[NSNumber numberWithInt:xxx]];
                } else {
                    NSMutableArray *newArray = [NSMutableArray arrayWithObject:[NSNumber numberWithInt:xxx]];
                    [dicF setObject:newArray forKey:kmer];
                }
                NSString *compKmer = complement(kmer);
                if ((array = [dicR objectForKey:compKmer])) {
                    // If complement kmer already exists, add position
                    [array addObject:[NSNumber numberWithInt:xxx + mer - 1]];
                } else {
                    NSMutableArray *newArray = [NSMutableArray arrayWithObject:[NSNumber numberWithInt:xxx + mer - 1]];
                    [dicR setObject:newArray forKey:compKmer];
                }
            }
        }
    }
    return dictionaries;
}

// Function to process TruSeq Standard mRNA Sample Prep kit data
// Read data is separated based on the strand it originates from
int main(int argc, const char * argv[]) {
    srand(time(NULL)); // Initialize with current time
    mer = 21;
    
    // Variables
    NSString *folderPathToStoreResult = @"path/to/yourFolder";
    NSString *outputFilePrefix = @"outputFilePrefix";
    NSArray *readPairs = [NSArray arrayWithObjects:
                            [NSArray arrayWithObjects:@"path/to/yourFile_R1.fastq",@"path/to/yourFile_R2.fastq", nil],
                          nil];
    NSString *referenceSeq = [NSString stringWithCString:"path/to/yourReferenceGenome.txt" encoding:NSUTF8StringEncoding];

    
    int totalPairCount = 0;
    int pairBothMappedCount = 0;
    int notValidPair = 0;
    int abnormalDistance = 0;
    int validDistancePair = 0;
    int totalRead = 0;

    // Read reference sequence
    NSString *fileContent = [NSString stringWithContentsOfFile:referenceSeq encoding:NSUTF8StringEncoding error:nil];
    // Hash reference sequence to create a dictionary of kmers and genomic locations
    headerAndSeq = [NSMutableDictionary dictionary];
    dictionayOfHeadderAndSeq(fileContent, headerAndSeq, 0);
    NSMutableDictionary *seqs_extendedIfNeeded = extendSequenceIfNeeded(headerAndSeq);
    dictionaries = getReadyForDictionary(seqs_extendedIfNeeded);
    NSMutableDictionary *superExtendSequence = extendSequence(headerAndSeq,2000);

    // Create, for each replicon, two int arrays for depth counting: one for the plus direction, one for the minus direction
    NSMutableDictionary *sNameAndPointersToCarray = [NSMutableDictionary dictionary];
    for (NSString *sName in seqs_extendedIfNeeded) {
        int size = [[headerAndSeq objectForKey:sName] length];
        int *pointerPlus = (int *)malloc(sizeof(int) * size);
        int *pointerMinus = (int *)malloc(sizeof(int) * size);
        for (int ccc = 0; ccc < size; ccc++) {
            pointerPlus[ccc] = 0;
            pointerMinus[ccc] = 0;
        }
        NSValue *p1 = [NSValue valueWithPointer:pointerPlus];
        NSValue *p2 = [NSValue valueWithPointer:pointerMinus];
        NSArray *array = [NSArray arrayWithObjects:p1,p2, nil];
        [sNameAndPointersToCarray setObject:array forKey:sName];
    }
    
    // Prepare array to store pair data objects
    NSMutableDictionary *repliconNamesAndPairData_all = [NSMutableDictionary dictionary];
    for (NSString *repliconName in dictionaries) {
        [repliconNamesAndPairData_all setObject:[NSMutableArray array] forKey:repliconName];
    }
    
    // Read and process pair data
    for (NSArray *array in readPairs) {
        NSString *f0 = array[0];
        NSString *f1 = array[1];
        FILE *fp1 = fopen([f0 cStringUsingEncoding:NSUTF8StringEncoding], "r");
        FILE *fp2 = fopen([f1 cStringUsingEncoding:NSUTF8StringEncoding], "r");
        int count = 0;
        char buf_read1_0[1024], buf_read1_1[1024], buf_read1_2[1024], buf_read1_3[1024];
        char buf_read2_0[1024], buf_read2_1[1024], buf_read2_2[1024], buf_read2_3[1024];
        while (1) {
            @autoreleasepool {
                if (count % 1000000 == 0) {
                    printf("%d\n",count);
                }
                if (count % 4 == 0) {
                    if (fgets(buf_read1_0, 512, fp1) == NULL) break;
                    if (fgets(buf_read2_0, 512, fp2) == NULL) break;
                } else if (count % 4 == 1) {
                    if (fgets(buf_read1_1, 512, fp1) == NULL) break;
                    if (fgets(buf_read2_1, 512, fp2) == NULL) break;
                } else if (count % 4 == 2) {
                    if (fgets(buf_read1_2, 512, fp1) == NULL) break;
                    if (fgets(buf_read2_2, 512, fp2) == NULL) break;
                } else if (count % 4 == 3) {
                    if (fgets(buf_read1_3, 512, fp1) == NULL) break;
                    if (fgets(buf_read2_3, 512, fp2) == NULL) break;
                }
                if (count % 4 == 3) {
                    totalPairCount++;
                    // Remove newline character from read1 sequence
                    buf_read1_1[strlen(buf_read1_1)-1] = '\0';
                    NSString *str1 = [NSString stringWithCString:buf_read1_1 encoding:NSUTF8StringEncoding];
                    // Remove newline character from read2 sequence
                    buf_read2_1[strlen(buf_read2_1)-1] = '\0';
                    NSString *str2 = [NSString stringWithCString:buf_read2_1 encoding:NSUTF8StringEncoding];
                    if ([str1 length] < mer || [str2 length] < mer) continue;
                    PairData *pd = validatePair2(str1, str2);
                    if (pd) {
                        if (pd.distance > 20 && pd.distance < 2000) {
                            NSMutableArray *xxx = [repliconNamesAndPairData_all objectForKey:pd.repliconName];
                            [xxx addObject:pd];
                            validDistancePair++;
                            NSArray *array = [sNameAndPointersToCarray objectForKey:pd.repliconName];
                            int *a;
                            if (pd.direction == 1) {
                                a = [[array objectAtIndex:0] pointerValue];
                            } else {
                                a = [[array objectAtIndex:1] pointerValue];
                            }
                            NSUInteger size = [[headerAndSeq objectForKey: pd.repliconName] length];
                            if (pd.direction == 1) {
                                if (pd.read1Location < pd.read2Location) {
                                    for (int i = pd.read1Location; i <= pd.read2Location; i++) {
                                        a[i]++;
                                    }
                                } else {
                                    for (int i = pd.read1Location; i < size; i++) {
                                        a[i]++;
                                    }
                                    for (int i = 0; i <= pd.read2Location; i++) {
                                        a[i]++;
                                    }
                                }
                            } else if (pd.direction == -1) {
                                if (pd.read1Location > pd.read2Location) {
                                    for (int i = pd.read2Location; i <= pd.read1Location; i++) {
                                        a[i]++;
                                    }
                                } else {
                                    for (int i = pd.read2Location; i < size; i++) {
                                        a[i]++;
                                    }
                                    for (int i = 0; i <= pd.read1Location; i++) {
                                        a[i]++;
                                    }
                                }
                            }
                        } else {
                            abnormalDistance++;
                        }
                    } else {
                        notValidPair++;
                    }
                }
                count++;
            }
        }
        fclose(fp1);
        fclose(fp2);
    }
    
    int thresholdDepthGT = 0;
    long sumOfReadLength = 0;
    
    // Output depth data for each replicon
    for (NSString *sname in sNameAndPointersToCarray) {
        NSArray *arr = [sNameAndPointersToCarray objectForKey:sname];
        int max = [[headerAndSeq objectForKey:sname] length];
        int *f = [[arr objectAtIndex:0] pointerValue];
        NSString *fileNamePlus = [NSString stringWithFormat:@"%@_%@_R.txt", outputFilePrefix, sname];
        NSString *filepathPlus = [folderPathToStoreResult stringByAppendingPathComponent:fileNamePlus];
        FILE *fp_Plus = fopen([filepathPlus cStringUsingEncoding:NSUTF8StringEncoding], "w");
        for (int mmm = 0; mmm < max; mmm++) {
            sumOfReadLength += f[mmm];
            if (f[mmm] > thresholdDepthGT) {
                NSString *data_Plus = [NSString stringWithFormat:@"%d\t%d\n", mmm+1, f[mmm]];
                fputs([data_Plus cStringUsingEncoding:NSUTF8StringEncoding], fp_Plus);
            }
        }
        fclose(fp_Plus);

        int *f2 = [[arr objectAtIndex:1] pointerValue];
        NSString *fileNameMinus = [NSString stringWithFormat:@"%@_%@_F.txt", outputFilePrefix, sname];
        NSString *filepathMinus = [folderPathToStoreResult stringByAppendingPathComponent:fileNameMinus];
        FILE *fp_Minus = fopen([filepathMinus cStringUsingEncoding:NSUTF8StringEncoding], "w");
        for (int mmm = 0; mmm < max; mmm++) {
            sumOfReadLength += f2[mmm];
            if (f2[mmm] > thresholdDepthGT) {
                NSString *data_Minus = [NSString stringWithFormat:@"%d\t%d\n", mmm+1, f2[mmm]];
                fputs([data_Minus cStringUsingEncoding:NSUTF8StringEncoding], fp_Minus);
            }
        }
        fclose(fp_Minus);
    }
    
    // Output depth and sequence data
    int windowSize = 120;
    int marginSize = 10;
    int stepSize = 10;
    for (NSString *sname in sNameAndPointersToCarray) {
        NSArray *arr = [sNameAndPointersToCarray objectForKey:sname];
        int repliconLength = [[headerAndSeq objectForKey:sname] length];
        int *f = [[arr objectAtIndex:0] pointerValue];
        int *f2 = [[arr objectAtIndex:1] pointerValue];
        NSString *fileNameFR = [NSString stringWithFormat:@"%@_%@_depthAndSeq3.txt", outputFilePrefix, sname];
        NSString *filepathFR = [folderPathToStoreResult stringByAppendingPathComponent:fileNameFR];
        FILE *fp_FR = fopen([filepathFR cStringUsingEncoding:NSUTF8StringEncoding], "w");
        
        for (int mmm = 0; mmm + windowSize < repliconLength; mmm += stepSize) {
            NSString *seq = [[superExtendSequence objectForKey:sname] substringWithRange:NSMakeRange(mmm,windowSize)];
            int sumOfDepth = 0;
            for (int kkk = marginSize; kkk < windowSize - marginSize; kkk++) {
                sumOfDepth += f[mmm+kkk];
                sumOfDepth += f2[mmm+kkk];
            }
            NSString *data_Plus = [NSString stringWithFormat:@"%d\t%@\n", sumOfDepth, seq];
            fputs([data_Plus cStringUsingEncoding:NSUTF8StringEncoding], fp_FR);
        }
        fclose(fp_FR);
    }
    
    // Output summary of analysis
    NSString *fileName5 = [NSString stringWithFormat:@"%@_summary.txt", outputFilePrefix];
    NSString *filepath5 = [folderPathToStoreResult stringByAppendingPathComponent:fileName5];
    NSMutableString *analysisReport = [NSMutableString stringWithFormat:@"Kmer: %d\n", mer];

    [analysisReport appendFormat:@"The depth greater than %d was exported.\n", thresholdDepthGT];
    [analysisReport appendFormat:@"Total number of read pairs:\t%d\n", totalPairCount];
    [analysisReport appendFormat:@"Number of pairs mapped in good distance:\t%d\n", validDistancePair];
    [analysisReport appendFormat:@"Number of pairs mapped in bad distance:\t%d\n", abnormalDistance];
    [analysisReport appendFormat:@"Number of pairs failed to be mapped:\t%d\n", notValidPair];
    [analysisReport appendFormat:@"Sum of read length mapped:\t%ld\n", sumOfReadLength];
    [analysisReport appendFormat:@"Read files loaded:\n"];
    for (NSArray *array in readPairs) {
        [analysisReport appendFormat:@"%@\t%@\n", array[0], array[1]];
    }
    [analysisReport appendFormat:@"References: \n"];
    for (NSString *seqName in headerAndSeq) {
        [analysisReport appendFormat:@"Sequence ID: %@ Length: %lu\n", seqName, (unsigned long)[[headerAndSeq objectForKey:seqName] length]];
    }
    [analysisReport writeToFile:filepath5 atomically:YES encoding:NSUTF8StringEncoding error:nil];
    
    return 0;
}

PairData *validatePair2(NSString *str1, NSString *str2) {
    @autoreleasepool {
        NSString *st1 = [str1 substringToIndex:mer];
        NSString *st2 = [str2 substringToIndex:mer];
        
        if ([st1 rangeOfString:@"N"].location != NSNotFound || [st2 rangeOfString:@"N"].location != NSNotFound) {
            // Return nil if sequences contain 'N'
            return nil;
        }

        NSMutableArray *allPossibilities = [NSMutableArray array];
        for (NSString *repliconName in dictionaries) {
            int len = [[headerAndSeq objectForKey:repliconName] length];
            NSArray *arr = [dictionaries objectForKey:repliconName];
            NSMutableDictionary *dicF = [arr objectAtIndex:0];
            NSMutableDictionary *dicR = [arr objectAtIndex:1];

            if ([dicF objectForKey:st1] && [dicR objectForKey:st2]) {
                NSMutableArray *arF = [dicF objectForKey:st1];
                NSMutableArray *arR = [dicR objectForKey:st2];
                for (NSNumber *numF in arF) {
                    int f = [numF intValue];
                    for (NSNumber *numR in arR) {
                        int r = [numR intValue];
                        PairData *pd = [[PairData alloc] init];
                        [allPossibilities addObject:pd];
                        pd.repliconName = repliconName;
                        pd.read1Location = f;
                        pd.read2Location = r;
                        
                        if (f < r) {
                            pd.distance = r - f + 1;
                        } else {
                            pd.distance = len - f + r;
                            pd.read2Location = r + len;
                        }
                        pd.direction = 1;
                    }
                }
            }

            if ([dicF objectForKey:st2] && [dicR objectForKey:st1]) {
                NSMutableArray *arF = [dicF objectForKey:st2];
                NSMutableArray *arR = [dicR objectForKey:st1];
                for (NSNumber *numF in arF) {
                    int f = [numF intValue];
                    for (NSNumber *numR in arR) {
                        int r = [numR intValue];
                        PairData *pd = [[PairData alloc] init];
                        [allPossibilities addObject:pd];
                        pd.repliconName = repliconName;
                        pd.read1Location = r;
                        pd.read2Location = f;
                        
                        if (f < r) {
                            pd.distance = r - f + 1;
                        } else {
                            pd.distance = len - f + r;
                            pd.read1Location = len + r;
                        }
                        pd.direction = -1;
                    }
                }
            }
        }

        if ([allPossibilities count] == 0) return nil;

        if ([allPossibilities count] == 1) {
            PairData *pd = [allPossibilities objectAtIndex:0];
            return pd;
        }

        NSArray *sortedArray = [allPossibilities sortedArrayUsingComparator:^(PairData *item1, PairData *item2) {
            int c1 = item1.distance;
            int c2 = item2.distance;
            return (NSComparisonResult)(c1 - c2);
        }];

        PairData *PairDataWithminmumDistance = [sortedArray objectAtIndex:0];
        int minimumDistance = PairDataWithminmumDistance.distance;
        NSMutableArray *arr = [NSMutableArray array];
        for (PairData *p in sortedArray) {
            if (p.distance == minimumDistance) {
                [arr addObject:p];
            }
        }
        int randamIndex = rand() % [arr count];
        return [sortedArray objectAtIndex:randamIndex];
    }
}

void dictionayOfHeadderAndSeq(NSString *string, NSMutableDictionary *dictionary, int mode) {
    // Convert mFASTA format data into dictionary with header as key and sequence as value, removing line breaks
    NSMutableString *mFASTAString = [NSMutableString stringWithString:string];
    [mFASTAString replaceOccurrencesOfString:@"\r" withString:@"\n" options:NSLiteralSearch range:NSMakeRange(0, [mFASTAString length])];
    NSArray *array = [mFASTAString componentsSeparatedByString:@">"];
    for (int aaa = 1 ; aaa < [array count]; aaa++) {
        NSString *oneFASTA = [array objectAtIndex:aaa];
        NSRange range = [oneFASTA rangeOfString:@"\n"];
        if (range.location == NSNotFound) continue;
        NSMutableString *header = [NSMutableString stringWithString:[oneFASTA substringToIndex:range.location]];
        if (mode == 0) {
            [header replaceOccurrencesOfString:@" " withString:@"_" options:NSLiteralSearch range:NSMakeRange(0, [header length])];
        } else if (mode == 1) {
            NSRange range2 = [header rangeOfString:@" "];
            if (range2.location != NSNotFound) {
                [header deleteCharactersInRange:NSMakeRange(range2.location,[header length] - range2.location)];
            }
        }
        NSMutableString *seq = [NSMutableString stringWithString:[oneFASTA substringFromIndex:range.location+1]];
        [seq replaceOccurrencesOfString:@"\n" withString:@"" options:NSLiteralSearch range:NSMakeRange(0, [seq length])];
        NSString *seqUpper = [seq uppercaseString];
        [dictionary setObject:seqUpper forKey:header];
    }
}

NSString *complement(NSString *DNA) {
    // Generate complementary DNA sequence
    const char *sequence = [DNA cStringUsingEncoding:NSUTF8StringEncoding];
    unsigned long length = strlen(sequence);
    char *complemet = (char *)malloc(sizeof(char) * length + 1);
    for (int ggg = 0 ; ggg < length ; ggg++) {
        complemet[ggg] = t_table[sequence[length - ggg - 1]];
    }
    complemet[length] = '\0';
    NSString *comp = [NSString stringWithCString:complemet encoding:NSUTF8StringEncoding];
    return comp;
}






