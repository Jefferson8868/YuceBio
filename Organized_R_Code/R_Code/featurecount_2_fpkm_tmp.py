def main():
    import pandas as pd
    import sys
    import argparse
    parser = argparse.ArgumentParser(
    '''
    -----------
    This is a little script to calculate the fpkm and tpm.
    The input file must be featurecounts result file
    -----------


    if you have any question please contact me,email:albert_xin@qq.com
    '''

    )
    parser.add_argument('-i','--file_in',
                       metavar = '\b', # this is a required para
                        help='please input your featurecounts result'
                       )
    parser.add_argument('-o','--file_out',
#                        meatvar = '\b' ,# this is a required para,
                        default = './data_out',
                        help = 'please set a out file, the default is data_out\fpkm\tpm'
                       )
    args = parser.parse_args()
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    else:
        pass

    file_in = args.file_in
    file_out = args.file_out
    fpkm_out = file_out+'fpkm.txt'
    tpm_out = file_out+'tpm.txt'


    # load featurecounts out data
    data_load = pd.read_csv(file_in,sep='\t',comment='#',index_col=0)

    # select the counts data
    data_felter = data_load.iloc[:,4:len(data_load)]

    # get the bamfile name
    data_list = list(data_felter.columns)[1:len(list(data_felter.columns))]

    # create new dataframe
    data_fpkm = pd.DataFrame(index=data_felter.index)

    # create  new dataframe
    data_tpm = pd.DataFrame(index=data_felter.index)


    # sum_length = sum(data_felter.Length)
    # calculate the fpkm
    for i in data_list:
#         print(i)
        sum_reads = sum(data_felter[i])
        data_felter['upper'] = data_felter[i]*10e9
        data_felter['downer']= sum_reads*data_felter.Length
        data_fpkm[i]=data_felter.upper/data_felter.downer

    #calculate the tpm
    for x in data_list:
        sum_fpkm = sum(data_fpkm[x])
        data_tpm[x] = (data_fpkm[x]*10e6)/sum_fpkm
    # export the data to csv
    data_fpkm.to_csv(fpkm_out,sep='\t')
    data_tpm.to_csv(tpm_out, sep='\t')



if __name__=='__main__':
    main()
