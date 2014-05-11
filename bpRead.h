#ifndef BPREAD_H
#define BPREAD_H
#include "main.h"
class bpString{
public:
	uint length;
	uint str_len;
	uchar * str;
	bpString(){
		str=NULL;
		length=0;
		str_len=0;
	}
	void initStr(uint src_length=0){
		if(str){
			delete [str_len]str;
		}
		length=src_length;
		str_len=(src_length+3)/4;
		if(str_len){
			str=new uchar[str_len];	
		}else{
			str=NULL;
		}
	}
	uchar getBP(uint index){
		uint str_idx=index/4;
		if(index>length || str_idx>str_len)throw ERROR_BPStringOutOfRange();		
		uchar bit_shift=(index%4)*2;
		uchar ch=(str[str_idx]>>(bit_shift))&0x03;
		return ch;
	}
	void writeBP(uint index,uchar bp){
		uint str_idx=index/4;
		uchar bit_shift=(index%4)*2;
		str[str_idx]=str[str_idx]&(~(0x03<<bit_shift))|(bp<<bit_shift);
	}
	void write(char * src,int len=-1,uint idx_begin=0){
		if(len==-1)len=strlen(src);
		initStr(len);
		for(uint i=0;i<length;i++){
			try{
				switch(src[i+idx_begin]){
				case 'A':writeBP(i,0);break;
				case 'T':writeBP(i,1);break;
				case 'C':writeBP(i,2);break;
				case 'G':writeBP(i,3);break;
				default:throw ERROR_BP_N();break;
				}
			}catch(ERROR_BP_N){
				cout<<"N";
				writeBP(i,0);
			}
		}
	}
	bpString(char * src){
		str=NULL;
		write(src);
	}
	bpString(char * src,uint len,uint idx_begin=0){
		str=NULL;
		write(src,len,idx_begin);
	}
	void copy(bpString * src){
		initStr(src->length);
		for(uint i=0;i<src->str_len;i++){
			str[i]=src->str[i];
		}
	}
	void print(){
		for(uint i=0;i<length;i++){
			switch(getBP(i)){
			case 0:cout<<'A';break;
			case 1:cout<<'T';break;
			case 2:cout<<'C';break;
			case 3:cout<<'G';break;
			default:throw ERROR_BP_N();break;
			}
		}
	}
	bool equalTo(bpString * target){
		if(target->length!=length)return false;
		for(int i=0;i<length;i++){
			if(target->getBP(i)!=getBP(i))return false;
		}
		return true;
	}
	/*
	int getOverlapLength(bString * target){
		int id_tail=length-1;
		int id_head=0;
		bool id_found=false;
		bool target_completed=false;
		while(id_tail>0){
			int id_tail_temp=id_tail;
			while((id_tail>=0) && (str[id_tail]!=target->str[id_head])){
				id_tail--;
			}
			if(id_tail<0)break; 
			id_tail_temp=id_tail;
			//记录暂时位点。从这里向后匹配。三种情况。
			//1.A串无法匹配到尾部，则放弃该匹配。
			//2.A串匹配到尾部。则接受该匹配。
			//3.B串被匹配至尾部，则认为B串是A串的子串。故应该丢弃B串。
			while(id_tail<length){
				if(str[id_tail]!=target->str[id_head])break;
				id_tail++;
				id_head++;
				if(id_head>=target->length){
					target_completed=true;
					break;
				}
				if(id_tail>=length){
					id_found=true;
					break;
				}
			}
			id_tail=id_tail_temp-1;
			id_head=0;
			if(id_found || target_completed){
				return length-id_tail_temp;
			}
		}
		return 0;
	}
	*/
};
class bpRead{
public:
	char seqID[READ_LEN_ID];
	bpString seq;
	char seq_src[READ_LEN_SEQ];
	char quality[READ_LEN_SEQ];
	uint pos_last;
	uint pos_cur;
	uchar status;
	void compileSeq(){
		seq.write(seq_src);
	}
	void preprocess(){

	}
	bpRead(){
		pos_last=0;
		pos_cur=0;
		status=STA_STOP;
	}
	uchar getBP(int idx){
		return seq.getBP(idx);
	}
	uchar getBPOffset(int offset){
		return seq.getBP(pos_cur+offset);
	}
	void pushForword(){
		pos_cur++;
		if(status==STA_START){
			pos_last=pos_cur;
		}else{
		}
		if(pos_cur>=seq.length-KMER_LENGTH+1){
			status=STA_END;
		}
	}
	bpString * getStr(int length){
		return new bpString(seq_src,length,pos_cur);
	}
	bool jointERROR(uint k=KMER_LENGTH){
		if(pos_cur-pos_last>k)return true;
		return false;
	}
	void start(){
		status=STA_START;
		pos_last=pos_cur;
	}
	void pause(){
		status=STA_PAUSE;
	}
	void stop(){
		status=STA_STOP;
		pos_cur=0;
		pos_last=0;
	}
};
#endif