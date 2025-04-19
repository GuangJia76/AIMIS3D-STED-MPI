import os, sys, inspect
from scipy.io.arff.arffread import r_relation
src_dir = os.path.dirname(inspect.getfile(inspect.currentframe()))
lib_dir = os.path.abspath(os.path.join(src_dir, 'lib'))
sys.path.insert(0, lib_dir)
import Leap
import socket
from queue import Queue
from iocp import Translate,Scale,Rotate

class SampleListener(Leap.Listener):
    
    def on_init(self, controller):
        print("Initialized")
        '''
        translate towards left or right
        '''
        self.tt_xPosition_pre = None
        self.left_right_queue_size = 64
        self.tt_left_right_queue = Queue(maxsize = self.left_right_queue_size)
        self.right = 0 # boolean 进 list 时会被转为 int，索性用 int
        self.l_index = 0
        self.left_right_queue_first = None
        '''
        translate towards front or backwards
        '''
        self.tt_zPosition_pre = None
        self.front_backwards_queue_size = 64
        self.tt_front_backwards_queue = Queue(maxsize = self.front_backwards_queue_size)
        self.front = 0 # boolean 进 list 时会被转为 int，索性用 int 
        self.f_index = 0
        self.front_backwards_queue_first = None
        '''
        translate towards up or down
        '''
        self.tt_yPosition_pre = None
        self.up_down_queue_size = 64
        self.tt_up_down_queue = Queue(maxsize = self.up_down_queue_size)
        self.up = 0 # boolean 进 list 时会被转为 int，索性用 int 
        self.u_index = 0
        self.up_down_queue_first = None
        '''
        scale 
        '''
        self.s_xPosition_pre = None
        self.scale_queue_size = 64
        self.scale_queue = Queue(maxsize = self.scale_queue_size)
        self.zoom_in = 0 # boolean  进list时会被转为int，索性用int 
        self.s_index = 0
        self.scale_queue_first = None
        '''
        rotate
        '''
        #left hand
        self.rot_left_tt_xPosition_pre = None
        self.rot_left_left_right_queue_size = 64
        self.rot_left_tt_left_right_queue = Queue(maxsize = self.rot_left_left_right_queue_size)
        self.rot_left_right = 0 # boolean 进 list 时会被转为 int，索性用 int
        self.rot_left_l_index = 0
        self.rot_left_left_right_queue_first = None
        
        self.rot_left_tt_zPosition_pre = None
        self.rot_left_front_backwards_queue_size = 64
        self.rot_left_tt_front_backwards_queue = Queue(maxsize = self.rot_left_front_backwards_queue_size)
        self.rot_left_front = 0 # boolean 进 list 时会被转为 int，索性用 int 
        self.rot_left_f_index = 0
        self.rot_left_front_backwards_queue_first = None
        
        self.rot_left_tt_yPosition_pre = None
        self.rot_left_up_down_queue_size = 64
        self.rot_left_tt_up_down_queue = Queue(maxsize = self.rot_left_up_down_queue_size)
        self.rot_left_up = 0 # boolean 进 list 时会被转为 int，索性用 int 
        self.rot_left_u_index = 0
        self.rot_left_up_down_queue_first = None
        
        #right hand
        self.rot_right_tt_xPosition_pre = None
        self.rot_right_left_right_queue_size = 64
        self.rot_right_tt_left_right_queue = Queue(maxsize = self.rot_right_left_right_queue_size)
        self.rot_right_right = 0 # boolean 进 list 时会被转为 int，索性用 int
        self.rot_right_l_index = 0
        self.rot_right_left_right_queue_first = None
        
        self.rot_right_tt_zPosition_pre = None
        self.rot_right_front_backwards_queue_size = 64
        self.rot_right_tt_front_backwards_queue = Queue(maxsize = self.rot_right_front_backwards_queue_size)
        self.rot_right_front = 0 # boolean 进 list 时会被转为 int，索性用 int 
        self.rot_right_f_index = 0
        self.rot_right_front_backwards_queue_first = None
        
        self.rot_right_tt_yPosition_pre = None
        self.rot_right_up_down_queue_size = 64
        self.rot_right_tt_up_down_queue = Queue(maxsize = self.rot_left_up_down_queue_size)
        self.rot_right_up = 0 # boolean 进 list 时会被转为 int，索性用 int 
        self.rot_right_u_index = 0
        self.rot_right_up_down_queue_first = None
        
        self.l_direction = None
        self.r_direction = None
        
        self.last_gesture = None #记录上一次的手势
        
    def on_connect(self, controller):
        print("Connected")
    
    def reset_rotate_direction(self):
        self.l_direction = None
        self.r_direction = None
        
    def set_translate_member_variable(self,list):
        '''
        translate towards left or right
        '''
        self.tt_xPosition_pre = list[0]
        self.left_right_queue_size = list[1]
        self.tt_left_right_queue = list[2]
        self.right = list[3]
        self.l_index = list[4]
        self.left_right_queue_first = list[5]
        '''
        translate towards front or backwards
        '''
        self.tt_zPosition_pre = list[6]
        self.front_backwards_queue_size = list[7]
        self.tt_front_backwards_queue = list[8]
        self.front = list[9]
        self.f_index = list[10]
        self.front_backwards_queue_first = list[11]
        '''
        translate towards up or down
        '''
        self.tt_yPosition_pre = list[12]
        self.up_down_queue_size = list[13]
        self.tt_up_down_queue = list[14]
        self.up = list[15] 
        self.u_index = list[16]
        self.up_down_queue_first = list[17]
        
    def set_rotate_left_hand_variable(self,list): 
        self.rot_left_tt_xPosition_pre = list[0]
        self.rot_left_left_right_queue_size = list[1]
        self.rot_left_tt_left_right_queue = list[2]
        self.rot_left_right = list[3] 
        self.rot_left_l_index = list[4]
        self.rot_left_left_right_queue_first = list[5]
        
        self.rot_left_tt_zPosition_pre = list[6]
        self.rot_left_front_backwards_queue_size = list[7]
        self.rot_left_tt_front_backwards_queue = list[8]
        self.rot_left_front = list[9] 
        self.rot_left_f_index = list[10]
        self.rot_left_front_backwards_queue_first = list[11]
        
        self.rot_left_tt_yPosition_pre = list[12]
        self.rot_left_up_down_queue_size = list[13]
        self.rot_left_tt_up_down_queue = list[14]
        self.rot_left_up = list[15] 
        self.rot_left_u_index = list[16]
        self.rot_left_up_down_queue_first = list[17]
        
    def set_rotate_right_hand_variable(self,list): 
        self.rot_right_tt_xPosition_pre = list[0]
        self.rot_right_left_right_queue_size = list[1]
        self.rot_right_tt_left_right_queue = list[2]
        self.rot_right_right = list[3]
        self.rot_right_l_index = list[4]
        self.rot_right_left_right_queue_first = list[5]
        
        self.rot_right_tt_zPosition_pre = list[6]
        self.rot_right_front_backwards_queue_size = list[7]
        self.rot_right_tt_front_backwards_queue = list[8]
        self.rot_right_front = list[9]
        self.rot_right_f_index = list[10]
        self.rot_right_front_backwards_queue_first = list[11]
        
        self.rot_right_tt_yPosition_pre = list[12]
        self.rot_right_up_down_queue_size = list[13]
        self.rot_right_tt_up_down_queue = list[14]
        self.rot_right_up = list[15]
        self.rot_right_u_index = list[16]
        self.rot_right_up_down_queue_first = list[17]        
        
    def on_disconnect(self, controller):
        print("Disconnected")

    def on_exit(self, controller):
        print("Exited")
    
    def send_socket(self,m):
        ip_port = ('127.0.0.1',6667)
        try:
            #请求连接服务端
            sk = socket.socket()
            sk.connect(ip_port)
            sk.send(m.encode())
            sk.close()
            print("发送成功")
        except Exception as e:
            print("发送失败")
            print(str(e))
    '''
    increasing 指代方向，当左右平移时，increasing=True指向右平移，其他同理；表示缩放时，increasing=True指放大
          数组中，下标在l_index之前的均满足递增或递减
    queue 首尾差距最大的为最终结果      
    '''
    def translate_or_scale(self,position,position_pre,position_queue_size,position_queue,increasing,l_index,queue_first):
        result=-1
        increasing = True if increasing==1 else False
        if Queue.qsize(position_queue) < position_queue_size:
            if(Queue.empty (position_queue)):
                position_pre=position
                queue_first=position  
            else:
                #if(not ((position_pre < position) ^  increasing)): #同或运算
                if(not (((position_pre < position) and increasing) or (((position_pre > position) and not increasing)))):                    
                    increasing = not increasing
                    while l_index>0:
                        position_queue.get()
                        l_index-=1    
                    queue_first=position  
            l_index+=1
            position_queue.put(position) 
        else:
            if(increasing):
                result=1
            else:
                result=0
                
            while Queue.qsize(position_queue) > 0:
                position_queue.get();
                l_index-=1
        increasing = 1 if increasing else 0
        return position_pre,position_queue_size,position_queue,increasing,l_index,result,queue_first
    '''
    get the most likely translate from "left、right、up、down、front、backwards" for given hand
    the parameter given is a list l_u_f_list which must follow the rule strictly
    '''
    def most_likely_translate(self,list,hand):
        abs_array=[-2.0,-2.0,-2.0] 
        
        l_r_position=hand.palm_position[0]
        list[0],list[1],list[2],list[3],list[4],l_r_result,list[5]=\
                                        self.translate_or_scale(l_r_position,list[0],\
                                        list[1],list[2],list[3],list[4],list[5])
        if not l_r_result == -1:
            abs_array[0]=abs(l_r_position-list[5])                   
       
        f_b_position=hand.palm_position[2]
        list[6],list[7],list[8],list[9],list[10],f_b_result,list[11]=\
                                        self.translate_or_scale(f_b_position,list[6],\
                                        list[7],list[8],list[9],list[10],list[11])
        if not f_b_result == -1:
            abs_array[1]=abs(f_b_position-list[11])
        
        u_d_position=hand.palm_position[1]
        list[12],list[13],list[14],list[15],list[16],u_d_result,list[17]=\
                                        self.translate_or_scale(u_d_position,list[12],\
                                        list[13],list[14],list[15],list[16],list[17])
        if not u_d_result == -1:
            abs_array[2]=abs(u_d_position-list[17])        
            
        if (l_r_result==-1 and f_b_result==-1 and u_d_result==-1) or max(abs_array)<30:            
            return False,None
        else:
            list[2].queue.clear()
            list[8].queue.clear()
            list[14].queue.clear()
            list[4]=0
            list[10]=0
            list[16]=0
            print("last:",self.last_gesture)
            index=abs_array.index(max(abs_array))
            if index==0:
                if l_r_result == 1:
                    return True,Translate.RIGHT                   
                elif l_r_result == 0:
                    return True,Translate.LEFT
            elif index==1:
                if f_b_result == 1:
                    return True,Translate.FRONT                   
                elif f_b_result == 0:
                    return True,Translate.BACKWARDS
            elif index==2:
                if u_d_result == 1 :
                    return True,Translate.UP
                elif u_d_result == 0 :                       
                    return True,Translate.DOWN
            return False,None
        
    def most_likely_rotate(self,left_hand,right_hand):
        l_list=[self.rot_left_tt_xPosition_pre,self.rot_left_left_right_queue_size,self.rot_left_tt_left_right_queue,self.rot_left_right,self.rot_left_l_index,self.rot_left_left_right_queue_first,\
              self.rot_left_tt_zPosition_pre,self.rot_left_front_backwards_queue_size,self.rot_left_tt_front_backwards_queue,self.rot_left_front,self.rot_left_f_index,self.rot_left_front_backwards_queue_first,\
              self.rot_left_tt_yPosition_pre,self.rot_left_up_down_queue_size,self.rot_left_tt_up_down_queue,self.rot_left_up,self.rot_left_u_index,self.rot_left_up_down_queue_first]
        l_result,l_direction=self.most_likely_translate(l_list,left_hand)        
        self.set_rotate_left_hand_variable(l_list)
        if l_result:
            self.l_direction = l_direction
        
        r_list=[self.rot_right_tt_xPosition_pre,self.rot_right_left_right_queue_size,self.rot_right_tt_left_right_queue,self.rot_right_right,self.rot_right_l_index,self.rot_right_left_right_queue_first,\
              self.rot_right_tt_zPosition_pre,self.rot_right_front_backwards_queue_size,self.rot_right_tt_front_backwards_queue,self.rot_right_front,self.rot_right_f_index,self.rot_right_front_backwards_queue_first,\
              self.rot_right_tt_yPosition_pre,self.rot_right_up_down_queue_size,self.rot_right_tt_up_down_queue,self.rot_right_up,self.rot_right_u_index,self.rot_right_up_down_queue_first]
        r_result,r_direction=self.most_likely_translate(r_list,right_hand)       
        self.set_rotate_right_hand_variable(r_list)
        if r_result:
            self.r_direction = r_direction
        
        if self.l_direction is not None and self.r_direction is not None:
            l_direction = self.l_direction
            r_direction = self.r_direction
            print("last:",self.last_gesture)
            
            if l_direction is Translate.UP and r_direction is Translate.DOWN:
                self.reset_rotate_direction()
                return True,Rotate.CLOCKWISE_Z
            elif l_direction is Translate.DOWN and r_direction is Translate.UP:
                self.reset_rotate_direction()
                return True,Rotate.COUNTERCLOCKWISE_Z
            elif l_direction is Translate.FRONT and r_direction is Translate.BACKWARDS:
                self.reset_rotate_direction()
                return True,Rotate.CLOCKWISE_Y
            elif l_direction is Translate.BACKWARDS and r_direction is Translate.FRONT:
                self.reset_rotate_direction()
                return True,Rotate.COUNTERCLOCKWISE_Y
            elif l_direction is Translate.LEFT and r_direction is Translate.RIGHT:
                self.reset_rotate_direction()
                return True,Scale.ZOOMIN
            elif l_direction is Translate.RIGHT and r_direction is Translate.LEFT:
                self.reset_rotate_direction()
                return True,Scale.ZOOMOUT
            
        return False,None
    
    
    def rotate_result(self,left_hand,right_hand):    
        result,action = self.most_likely_rotate(left_hand,right_hand)
        if result:
            #去除做完动作后的返回动作 
            if (action is Rotate.CLOCKWISE_Z        and self.last_gesture is Rotate.COUNTERCLOCKWISE_Z ) or\
               (action is Rotate.COUNTERCLOCKWISE_Z and self.last_gesture is Rotate.CLOCKWISE_Z        ) or\
               (action is Rotate.CLOCKWISE_X        and self.last_gesture is Rotate.COUNTERCLOCKWISE_X ) or\
               (action is Rotate.COUNTERCLOCKWISE_X and self.last_gesture is Rotate.CLOCKWISE_X        ) or\
               (action is Rotate.CLOCKWISE_Y        and self.last_gesture is Rotate.COUNTERCLOCKWISE_Y ) or\
               (action is Rotate.COUNTERCLOCKWISE_Y and self.last_gesture is Rotate.CLOCKWISE_Y        ):
                self.last_gesture = None
            else:
                self.last_gesture = action
            if  self.last_gesture is not None:    
                self.send_socket(str(action.value))
                print(action)
        
    def scale_result(self,s_c_position):#以两手之间的距离差值进行判断
        self.s_xPosition_pre,self.scale_queue_size,self.scale_queue,self.zoom_in,self.s_index,s_result,self.scale_queue_first=\
                                         self.translate_or_scale(s_c_position,self.s_xPosition_pre,\
                                        self.scale_queue_size,self.scale_queue,self.zoom_in,self.s_index,self.scale_queue_first)
        if s_result == 1:
            print("放大")
            #self.send_socket(str(Scale.ZOOMIN.value))
        elif s_result == 0:
            print("缩小")
            #self.send_socket(str(Scale.ZOOMOUT.value))
        
    def translate_result(self,hand):
        list=[self.tt_xPosition_pre,self.left_right_queue_size,self.tt_left_right_queue,self.right,self.l_index,self.left_right_queue_first,\
              self.tt_zPosition_pre,self.front_backwards_queue_size,self.tt_front_backwards_queue,self.front,self.f_index,self.front_backwards_queue_first,\
              self.tt_yPosition_pre,self.up_down_queue_size,self.tt_up_down_queue,self.up,self.u_index,self.up_down_queue_first]
        result,direction=self.most_likely_translate(list,hand)
        self.set_translate_member_variable(list)
        
        if result:
            #去除做完动作后的返回动作 
            if (direction is Translate.RIGHT      and self.last_gesture is Translate.LEFT      ) or \
               (direction is Translate.LEFT       and self.last_gesture is Translate.RIGHT     ) or \
               (direction is Translate.FRONT      and self.last_gesture is Translate.BACKWARDS ) or \
               (direction is Translate.BACKWARDS  and self.last_gesture is Translate.FRONT     ) or \
               (direction is Translate.UP         and self.last_gesture is Translate.DOWN      ) or \
               (direction is Translate.DOWN       and self.last_gesture is Translate.UP        ):
                self.last_gesture = None
            else:
                self.last_gesture = direction
            if  self.last_gesture is not None:
                self.send_socket(str(direction.value))
                print(direction)
            
    def on_frame(self, controller):
       
        frame = controller.frame()
        try:
            #print(frame)
            if len(frame.hands)==2:
                left_hand = frame.hands[0] if frame.hands[0].is_left else frame.hands[1]
                right_hand = frame.hands[0] if frame.hands[0].is_right else frame.hands[1]
                self.rotate_result(left_hand,right_hand)
            elif len(frame.hands)==1:
                pass
                self.translate_result(frame.hands[0])
                
        except Exception as e:
            print(e.message)
        

def main():    
    listener = SampleListener()
    controller = Leap.Controller()
    controller.add_listener(listener)
    
    print("input exit to quit...")
    try:
        con=sys.stdin.readline()
        while not con=='exit':
            con=sys.stdin.readline()
    except KeyboardInterrupt:
        pass
    finally:
        controller.remove_listener(listener)
    #socket断开时退出本程序    
if __name__ == "__main__":
    main()


