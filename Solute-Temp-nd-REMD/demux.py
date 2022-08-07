import copy


f = open('rem.log','r')
started =0
out=[0,1,2,3,4, 1]
prev=[0,1,2,3,4 ,1 ]
outfile= open ('outdemux.dat' ,'w')

outfile.writelines( '#exchange, traj1, traj2, traj3, traj4, position of ham0')

for line in f.readlines()[1:100]:
    #print(started , prev ,out, int(line.split()[0]) , out[int(line.split()[0])], prev[0])
    if line[:10] != '# exchange':
        if started ==0 :
            continue
        else:
            if line.split()[7]== 'F':

                out[int(line.split()[0])]=prev[int(line.split()[0])]
            else :
                out[int(line.split()[0])]=prev[int(line.split()[1])]
    else  :
            if started ==1 :
                out[5]=out.index(1)
                outfile.writelines( '%s %s %s %s %s                %s \n'  %(out[0] ,out[1], out[2] , out[3] , out[4] , out[5]))
                prev = copy.deepcopy(out)

                out[0] =f'{int(line.split()[2]):05}'#print("{:02d}".format(1))
            else :
                started =1
                out[0] = '00000'
