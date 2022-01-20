classdef neuron
   properties
        a ;
        b ;
        c ;
        d ;
    
   end
   
   methods
       
       function obj = neuron(a,b,c,d)
           obj.a = a;
           obj.b = b;
           obj.c = c;
           obj.d = d;
            
       end
       function obj = change_d(obj,d_value)
           
            obj.d = d_value;
           
       end
    end
end