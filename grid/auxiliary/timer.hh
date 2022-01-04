#ifndef TIMER_HH
#define TIMER_HH

template< class Decomposition >

///@cond

class Timer {
public:
    Timer( Decomposition & decomp ) : decomp_(decomp) {}
    
    void begin() { start_time = clock(); }
    
    void stop_local(const char * str) {
        double time_local = (float) (clock() - start_time) / CLOCKS_PER_SEC;
        if( decomp_.master() )
            std::cout << "time for: " << str << " || local: " << time_local << std::endl; 
    }
    
    void stop(const char * str) {
        double time_local = (float) (clock() - start_time) / CLOCKS_PER_SEC;
        double time_global = time_local;
        double time_max = time_local;
        decomp_.scalar_sum( time_global );
        decomp_.scalar_max( time_max );
        if( decomp_.master() )
            std::cout << "time ( local || average || max ): " << str << "\t \t \t || \t\t" << time_local << "\t" << time_global / decomp_.size()  << "\t" << time_max <<  std::endl;
    }
    
    void getTime(double &time_local, double &time_average, double &time_max) {
        time_local = (float) (clock() - start_time) / CLOCKS_PER_SEC;
        double time_global = time_local;
        time_max = time_local;
        decomp_.scalar_sum( time_global );
        decomp_.scalar_max( time_max );
        time_average = time_global / decomp_.size();
    }
    
    void getTime_local( double &time_local ) {
        time_local = (float) (clock() - start_time) / CLOCKS_PER_SEC;
    }
    
private:
    Decomposition &decomp_;
    double start_time;
};

class TimerLocal {
public:
    void start() { begin(); }
    void begin() { start_time = clock(); }
    
    void stop(const char * str) {
        double time_local = (float) (clock() - start_time) / CLOCKS_PER_SEC;
        std::cout << "time for: " << str << " || local: " << time_local << std::endl; 
    }
            
private:
    double start_time;
};

class TimerEmpty {
public:
    void begin() {}    
    void stop(const char *) {}
};

///@endcond

#endif // TIMER_HH