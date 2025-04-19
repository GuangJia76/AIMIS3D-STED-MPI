数据库安装指南
1.安装mysql
2.安装mysql_odbc连接器
3.电脑 控制面板->管理工具->odbc数据源配置 （详情可见百度：关键字 如何配置数据库ODBC数据源）
4.代码db_load.pro 里 的INITDatabase 过程中，更改下面的三个字符串变量的值为数据源中对应的值即可。 
	datasourse
	use_id
	password

如果INITDatabase 函数运行成功，则表示数据库连接成功。

5.在mysql 命令台或者管理界面打开 database 文件夹下的create.sql 文件，并运行，运行成功则表示数据库创建成功。

如有疑问欢迎来信询问
邮箱：632446790@qq.com
 
 