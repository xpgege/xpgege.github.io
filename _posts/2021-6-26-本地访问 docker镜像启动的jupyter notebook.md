---
layout:     post
title:      本地访问 docker镜像启动的jupyter notebook
subtitle:   docker/jupyter notebook
date:       2021-6-26
author:     xpgege
header-img: img/post-web.jpg
catalog: true
tags:
   - 知识记录
---

#1.登录服务器，安装好jupyter notebook，生成配置文件：

```
$ jupyter notebook --generate-config
Writing default config to: /home/xxx/.jupyter/jupyter_notebook_config.py
# jupyter_notebook_config.py 就是配置文件
```

#2.使用 Ipython 生成密码
```
In [1]: from notebook.auth import passwd

In [2]: passwd()
Enter password:
Verify password:
Out[2]: 'sha1:*********************************c'
# Enter password 是在本地浏览器登陆使用的密码，Out[2] 输出的是填写在配置文件中的密钥

In [3]: exit()
```

#3.编辑配置文件
**注意配置文件全部都是注释掉的，记得取消注释呀，555被坑了好久**
![image.png](https://i.loli.net/2021/06/26/EO6aQu9HTqFkhyd.png)
```
# 在该文件中搜索以下配置，没有的话就添加


c.NotebookApp.allow_remote_access = True    # 可能没有需要自己添加
c.NotebookApp.ip = '*'
c.NotebookApp.open_browser = False
c.NotebookApp.password = u'sha1:*********************************c' # 填写刚才在第二步Out[2]生成的密钥
```

#4.运行
在服务器运行 jupyter notebook。 然后在本地打开浏览器，网址处输入hostname:port，其中hostname就是服务器的ip地址，例如，222.20.79.247:8888，然后输入在第二步设置的密码即可。
![image.png](https://i.loli.net/2021/06/26/f1Awt9XR3VmOsqa.png)

**踩坑+1**
由于我是从docker镜像启动jupyter notebook
```
sudo  docker run --rm -v $(pwd):/data -p 8888:8888 --entrypoint="" kavrakilab/hla-arena jupyter notebook --port=8888  --ip=0.0.0.0 --allow-root
```
此时我本地浏览器登陆222.20.79.247:8888
登陆界面为
![image.png](https://i.loli.net/2021/06/26/LeVJYUDXE4QK2vR.png)

此时我输入当初设置的密码为invalid，而token令牌开始我也不知道
后来发现每次登陆jupyter 都会随机产生一个token，将以下token输入即课正常登陆
![image.png](https://i.loli.net/2021/06/26/va2Npd8YkoUgKMC.png)