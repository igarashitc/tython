## memorondom for matplotlib

#基本的なこと
fig = plt.figure()
ax  = plt.subplot()
#or 
fig, ax = plt.subplots()

#軸のあたいを x10^*のような書式にする
ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax.ticklabel_format(style='sci', axis='x', scilimit=(0,0))

#colorbarの設定（仮）
#pcolormeshとかにcolorbarをつける
pcl = ax.pcolormesh(x,z,data)
ax.colorbar(pcl, ax=ax)



