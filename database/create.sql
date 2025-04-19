#drop database if exists test;
create database if not exists test CHARACTER SET UTF8;
#create the table;
use test;

drop table if exists data_table;
drop table if exists case_info;
drop table if exists patient_info;
create table IF NOT exists patient_info
(id int auto_increment primary key,
patient_name varchar(10),
gender enum('male','female') not null default 'male' ,
birth date,
age int unsigned,
patient_ID_card char(18) unique
);

CREATE TABLE IF NOT EXISTS case_info (
    id INT AUTO_INCREMENT PRIMARY KEY,
    modaling ENUM('CT', 'MR', 'DSA', 'MRI'),
    hospital VARCHAR(20),
    exam_date DATETIME,
    doctor_name VARCHAR(10),
    description VARCHAR(50),
    note VARCHAR(50),
    path VARCHAR(300),
    _snapshot VARBINARY(8),
    _data VARBINARY(10),
    patient_id INT,
    FOREIGN KEY (patient_id)
	REFERENCES patient_info (id)
);

CREATE  TABLE IF NOT exists data_table(
		id INT PRIMARY KEY,
		_snapshot BLOB(104857600),
		imr int,
        imc int,
		_data BLOB(104857600),
        dr int,
        dc int,
        dd int,
        FOREIGN KEY (id)
        REFERENCES case_info(id)
);



drop procedure if exists insert_info;
DELIMITER $
#create procedure
create procedure  insert_info(IN patient_name varchar(10),
IN str_gender varchar(10),
IN str_birth varchar(10),
IN age int unsigned,
IN patient_ID_card char(18),
IN str_modaling varchar(10),
IN hospital VARCHAR(20),
IN exam_date DATETIME,
IN doctor_name VARCHAR(10),
IN description VARCHAR(50),
IN note VARCHAR(50),
IN path VARCHAR(300)
)
BEGIN
		declare patient_id int;
        declare formatstr varchar(10);
		declare birth datetime;
        declare gender enum('male','female');
        declare modaling ENUM('CT', 'MR', 'DSA','MRI');
        set formatstr='%Y%m%d';
        set birth=str_to_date(str_birth,formatstr);
        set gender= case when str_gender='' then 'male' else str_gender end;
        set modaling = case when str_modaling ='' then null else str_modaling end;
        
        #case : no id card 
		if not exists(select * from patient_info as a where 
			a.patient_name = patient_name and
			(case when isnull(gender) then 1=1 else	a.gender=gender end) and
		   (case when isnull(birth) then 1=1 else	a.birth=birth end) and
		   (case when isnull(age) then 1=1 else	a.age=age end) and
		   (case when isnull(patient_ID_card) then 1=1 else	 a.patient_ID_card=patient_ID_card end)
            )then
            
			insert into patient_info(id,patient_name,gender,birth,age,patient_ID_card) values(null
            ,patient_name
            ,gender
            ,birth
            ,case when isnull(age) and str_birth !='' then (year(now())-year(birth)-1)+ (date_format(birth,'%m%d')<=DATE_FORMAT(NOW(),'%m%d')) else age end
            ,patient_ID_card); 
		end if;
		set patient_id =(select max(id)  from patient_info as a where 
			a.patient_name = patient_name and
           (case when isnull(gender) then 1=1 else	a.gender=gender end) and
		   (case when isnull(birth) then 1=1 else	a.birth=birth end) and
		   (case when isnull(age) then 1=1 else	a.age=age end) and
		   (case when isnull(patient_ID_card) then 1=1 else	 a.patient_ID_card=patient_ID_card end)
           );
       #insert the case 
		insert into case_info(id,modaling,hospital,exam_date,doctor_name,description,note,path,patient_id) 
		values(null,modaling,hospital,exam_date,doctor_name,description,note,path,patient_id);
        
       # select id,exam_date from case_info where patient_id =patient_id  order by id desc limit 1;
       
       select id from case_info as a where
		a.exam_date=exam_date and
        a.patient_id=patient_id and
		(case when isnull(modaling) then 1=1 else a.modaling= modaling end) and
        (case when isnull(hospital) then 1=1 else a.hospital= hospital end) and
        (case when isnull(doctor_name) then 1=1 else a.doctor_name=doctor_name end) and
        (case when isnull(description) then 1=1 else a.description=description end) and
        (case when isnull(path) then 1=1 else a.path=path end) and
        (case when isnull(note) then 1=1 else a.note=note end)
        order by id desc limit 1;
       
       
END
$
DELIMITER ;







/*
drop procedure if exists select_data;
DELIMITER $
#create procedure
create procedure  select_data(in id int)
BEGIN
		insert into  data_table(_snapshot,_data) select _snapshot,_data from  case_info where  case_info.id=id;
END
$
DELIMITER ;



*/



drop view if exists show_view;
create view show_view as
select b.id as id 
,patient_name as name
,gender
,age
,birth
,modaling
,exam_date as date
,hospital
,doctor_name as doctor
,description
,note
,path
,c._snapshot as image
,c._data as biData
from patient_info  as a inner join case_info as b on a.id=b.patient_id left join data_table as c on b.id=c.id ;


drop procedure if exists search_info;
DELIMITER $
#create procedure
create procedure  search_info(IN pname varchar(10),
IN modaling varchar(10),
IN date1 varchar(10),
IN date2 varchar(10)
)
BEGIN
		declare formatstr varchar(10);
		declare d1 datetime;
        declare d2 datetime;
        
        set formatstr='%Y%m%d';
        set d1=str_to_date(date1,formatstr);
        set d2=str_to_date(date2,formatstr);
        
        select id,name,gender,age,date_format(birth,'%Y-%m') as birth,a.modaling,date_format(date,'%Y-%m-%d %h:%i:%s') as date,hospital,doctor from show_view as a where
        (case when pname='' then 1=1 else a.name=pname end) and
        (case when modaling = '' then 1=1 else a.modaling =modaling end) and
        (case when date1='' and date2='' then 1=1  
			  when date1='' then a.date between d2 and DATE_ADD(d2,INTERVAL 1 DAY)
			  when date2='' then a.date between d1 and DATE_ADD(d1,INTERVAL 1 DAY)
		else
			  a.date between d1 and d2 
		end
		);
		
END
$
DELIMITER ;

drop procedure if exists delete_info;
DELIMITER $
#create procedure
create procedure  delete_info(in _id int)
BEGIN
	 declare  pid int;
	 set pid = (select patient_id from case_info where id=_id);
     delete from data_table where id=_id;
     delete from case_info  where id=_id;
     if not exists(select * from case_info where patient_id =pid) then 
		delete from patient_info where id=pid;
     end if; 
		
END
$
DELIMITER ;


drop procedure if exists update_info;
DELIMITER $
#create procedure
create procedure  update_info(IN case_id int,
IN patient_name varchar(10),
IN str_gender varchar(10),
IN str_birth varchar(10),
IN patient_ID_card char(18),
IN str_modaling varchar(10),
IN hospital VARCHAR(20),
IN doctor_name VARCHAR(10),
IN description VARCHAR(50),
IN note VARCHAR(50),
IN path VARCHAR(300)
)
BEGIN
		declare formatstr varchar(10);
        set formatstr='%Y%m%d';
		update patient_info as a
        set patient_name=patient_name
        ,birth=str_to_date(str_birth,formatstr)
        ,age=case when str_birth !='' then (year(now())-year(birth)-1)+ (date_format(birth,'%m%d')<=DATE_FORMAT(NOW(),'%m%d')) else 0 end
		,gender= case when str_gender='' then 'male' else str_gender end 
		where a.id =(select b.patient_id from case_info as b where b.id=case_id limit 1);

        update case_info
        set modaling = case when str_modaling ='' then null else str_modaling end
        ,hospital =hospital
		,doctor_name= doctor_name
        ,description=description
        ,note=note
        ,path=path
        where id=case_id;
END
$
DELIMITER ;



drop procedure if exists data_search;
DELIMITER $
#create procedure
create procedure  data_search(IN case_id int)
BEGIN
	select _snapshot as image,imr as imRow,imc as imCol,_data as data ,dr as dataRow,dc as dataCol,dd as dataDim from data_table where id=case_id;
END
$
DELIMITER ;



drop procedure if exists id_search;
DELIMITER $
#create procedure
create procedure  id_search(IN case_id int)
BEGIN
	select name,gender,date_format(birth,'%Y') as y,date_format(birth,'%m') as m,date_format(birth,'%d') as d,modaling,hospital,doctor,description,note,path from show_view as a where id=case_id; 
END
$
DELIMITER ;


