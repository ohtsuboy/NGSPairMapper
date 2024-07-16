//
//  PairData.h
//  analyzeTag1
//
//  Created by 大坪 嘉行 on 2016/05/08.
//  Copyright © 2016年 大坪 嘉行. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface PairData : NSObject


@property NSString *repliconName;
@property NSString *read1_genomeSeq;
@property NSString *read2_genomeSeq;
@property int read1Location;
@property int read2Location;
@property int distance;
@property int direction;





@end
