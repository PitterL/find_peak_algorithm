// ConsoleFindpeaks.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <iostream>
#include <vector>

using namespace std;

struct peak{
	int val;	//peak value
	int index;	//peak channel number
	int sign;	//second derivatvie value
};

typedef unsigned char u8;
struct position{
	u8 index;	//channel number integer
	u8 index_decimal;	//chnnel number decimal
	int val;	//value
};

enum{
	START,
	END,
	RANGE
};

struct hamp{
	vector<struct peak>::iterator full[RANGE];	//store hamp channel range
	struct position point;
};

struct merge_config{
	int thld;	//condition 2: positive peak value must less than threshold 
	int hysis;	//condition 1: to adjcent high peak less than hysis
	int count;	//condition 1: OR two positie peak distance less than count
};

//diff will ignore hysis ripple
int inline caculate_diff(int val, int prev, int hysis)
{
	int diff;

	diff = val - prev;
	if (abs(diff) < hysis)
		diff = 0;

	return diff;
}

void findPeaks(int *buf, int count, int hysis, vector<struct peak> &vector_peak)
{
	/*	i: the search squerence number
		j: offset to caculate buf[i] - buf[i - 1]
		k: offset to caculate buf[i + 1] - buf[i]
	*/
	int i, j, k;
	/*	diff[0]: buf[i] - buf[i - 1]
		diff[1]: buf[i + 1] - buf[i]
	*/
	int diff[2];
	struct peak val_p;

	for (i = 1, j = 0, k = 1; i + k < count;)	{

		//WARNING_ON(j >= k);

		//calulate fist derivative
		diff[0] = caculate_diff(buf[i + j], buf[i - 1], hysis);
		diff[1] = caculate_diff(buf[i + k], buf[i + j], hysis);

		//printf("[%d %d %d (%d %d %d)] = %d %d %d [%d %d]\n", i - 1, i + j, i + k, i, j, k, buf[i - 1], buf[i + j], buf[i + k], diff[0], diff[1]);

		if (diff[0] && diff[1]) {
			if ((diff[0] ^ diff[1]) < 0) {   //first derivative opposite sign, there will be a maximum/minimum
				val_p.index = i + j /*+ ((k-j) >> 1)*/;	//if index use i + j: first index of maximum/minimus in same value
													//if index use j + j + ((k-j)>>1): middle index of maximum/minimus in same value
				val_p.val = buf[val_p.index];
				val_p.sign = diff[1] - diff[0];  //second derivative
				vector_peak.push_back(val_p);

				//printf("*%d: %d %s\n", val_p.index, val_p.val, val_p.sign > 0 ? "MIN":"MAX");
			}

			//move index after i + j
			i += j + 1;	//combine j to i, so k will be re-caculated for i + k offset same
			k -= j;
			j = 0;
		}
		else {
			if (!diff[0]) {		//re-search valid diff[0]
				for (j = 1 ; i + j < count - 1; j++) {		//search j, so end before count - 1
					diff[0] = caculate_diff(buf[i + j], buf[i - 1], hysis);
					if (diff[0])
						break;
				}
			}

			if (diff[0]) {		//re-search valid diff[1]
				for (k = j + 1; i + k < count; k++) {
					diff[1] = caculate_diff(buf[i + k], buf[i + j], hysis);
					if (diff[1])
						break;
				}
			}

			if (!diff[0] || !diff[1])	//no valid devriative could be found
				break;
		}
	}
}

void append_peaks_in_edge(int *buf, int count, vector<struct peak> &vector_peak, int thld)
{
	vector<struct peak>::iterator iter;
	struct peak val_p;

	//add peak in front
	iter = vector_peak.begin();
	if (iter->index != 0) {
		val_p.index = 0;
		val_p.val = buf[val_p.index];
		if (iter->sign > 0 &&
			iter->val + thld <= buf[0]) {	//Maximum
			val_p.sign = -1;
		}
		else if (iter->sign < 0 &&
			iter->val - thld >= buf[0]) {	//Minimum
			val_p.sign = 1;
		}
		else {
			val_p.sign = 0;
		}
		vector_peak.insert(vector_peak.begin(), val_p);
	}

	//add peak in end
	iter = vector_peak.end() - 1;
	if (iter->index != count - 1) {
		val_p.index = count - 1;
		val_p.val = buf[val_p.index];
		
		if (iter->sign > 0 &&
			iter->val + thld <= buf[count - 1]) {	//Maximum
			val_p.sign = -1;
		}
		else if (iter->sign < 0 &&
			iter->val - thld >= buf[count - 1]) {	//Minimum
			val_p.sign = 1;
		}
		else{
			val_p.sign = 0;
		}

		vector_peak.push_back(val_p);
	}
}

void merge_adjacent_high_peak(vector<struct peak> &vector_peak, int thld, const struct merge_config &mcfg, vector<struct hamp> &vector_hamp)
{
	/* iter:
	0 : first high peak
	1 : second high peak
	2 : low peak after iter[0]
	3 : low peak before iter[1], --- may equal iter[2]*/
	vector<struct peak>::iterator iter0, iter1, iter2, iter3, iter_next, iter_prev, iter_curr;
	int distance,diff;
	struct hamp h_val;

	for (iter0 = vector_peak.begin(); iter0 != vector_peak.end();) {

		if (iter0->sign < 0 && iter0->val > thld) {		//check maximum threshold
			if (iter0 > vector_peak.begin())
				iter_prev = iter0 - 1;
			else
				iter_prev = iter0;
			h_val.full[START] = iter_prev;

			if (iter0 < (vector_peak.end() - 1))
				iter_next = iter0 + 1;
			else
				iter_next = iter0;
			h_val.full[END] = iter_next;

			for (iter1 = iter0 + 1; iter1 != vector_peak.end(); iter1++) {
				if (iter1->sign < 0 && iter1->val > thld) {		//check maximum threshold
					distance = iter1->index - iter0->index;	//it's a positive value
					diff = iter0->val - iter1->val;
					diff = abs(diff);
					if (distance < mcfg.count || diff < mcfg.hysis) {	//merge mandotory condition 1: distance < merge count OR diff < hysis
						iter2 = iter0 + 1;
						iter3 = iter1 - 1;

						if (iter2->sign > 0 && iter3->sign > 0) {	//these two value must > 0, it's minimum
							if (iter0->val - iter2->val < mcfg.thld ||
								iter1->val - iter3->val < mcfg.thld) {	//merge mandotory condition 2: peak thld < merge threshold
								if (iter0 < (vector_peak.end() - 1))
									iter_next = iter1 + 1;
								else
									iter_next = iter1;
								h_val.full[END] = iter_next;
							}
							else
								break;
						}
					}
					else
						break;
				}
			}

			iter0 = iter1;
			vector_hamp.push_back(h_val);
		}
		else
			iter0++;
	}
}

//range find: from first data large than threshold to last data large than threshold
int cacluate_range(const int *buf, int st, int len, int thld)
{
	int i, j;

	for (i = 0; i < len; i++) {
		if (buf[st + i] > thld) {
			for (j = len - 1; j >= i; j--) {
				if (buf[st + j] > thld) {
					return (((st + i) << 16) | (j - i + 1));
				}
			}
		}
	}

	return 0;
}

//average value
int caculate_average_value(const int *buf, int st, int len)
{
	int i, sum, avg;

	sum = avg = 0;
	for (i = 0; i < len; i++) {
		sum += buf[st + i];
	}

	if (sum) {
		avg = sum / len;
	}

	return avg;
}

//variance value
int caculate_variance_value(const int *buf, int st, int len)
{
	int i, avg, e_1, e_sum;

	avg = caculate_average_value(buf, st, len);

	e_sum = 0;
	for (i = 0; i < len; i++) {
		e_1 = buf[st + i] - avg;
		e_1 *= e_1;
		e_sum += e_1;
	}

	return e_sum;
}

//position: from middle, check strength left and right, caculate percent 
void caculate_position(const int *buf, int st, int len, struct position &pos)
{
	int i, mid, sum[2], percent;

	pos.index = 0;
	pos.index_decimal = 0;
	pos.val = 0;

	if (len == 1) {
		pos.val = buf[st];
		pos.index = st;
	}else {
		sum[0] = sum[1] = 0;
		mid = len >> 1;
		for (i = 0; i < mid; i++) {
			sum[0] += buf[st + i];
			sum[1] += buf[st + len - i - 1];
		}
		if (sum[0] + sum[1]) {
			percent = (sum[0] * 100) / (sum[0] + sum[1]);
			pos.val = (buf[st + mid - 1] * percent + buf[st + ((len + 1) >> 1)] * (100 - percent)) / 100;
			pos.index = st + mid - 1;
			pos.index_decimal = 100 - percent;
		}
	}
}

//search where is the middle point of each hamp
/*  1 caculate range which data value higher than average value
	2 higher average value
	3 caculate again
	exit condition: <1> array len less than 3 <2> average steady <3>  variance less than hysis */
void caculate_peak_in_hamp(int *buf, vector<struct hamp> &vector_hamp, int thld, int hysis)
{
	vector<struct hamp>::iterator h_iter;
	vector<struct peak>::iterator p_iter;
	struct position pos;
	int val, st[2], len[2], avg[2], avg_diff, e_sum, e_diff;

	for (h_iter = vector_hamp.begin(); h_iter != vector_hamp.end(); h_iter++) {

		st[0] = st[1] = h_iter->full[START]->index;
		len[0] = len[1] = h_iter->full[END]->index - st[1] + 1;
		avg[0] = avg[1] = thld;

		do {
			val = cacluate_range(buf, st[1], len[1], avg[1]);
			if (val) {
				avg[0] = avg[1];
				st[0] = st[1];
				len[0] = len[1];

				st[1] = (val >> 16);
				len[1] = val & 0xffff;
				avg[1] = caculate_average_value(buf, st[1], len[1]);
				e_sum = caculate_variance_value(buf, st[1], len[1]);

				avg_diff = avg[0] - avg[1];
				e_diff = hysis * hysis;

				//printf("%d %d: %d(%d, %d) %d\n", st[1], len[1], avg_diff, avg[0], avg[1], e_sum);
				if (len[1] <= 3 || e_sum <= e_diff) {
					st[0] = st[1];
					len[0] = len[1];
					break;
				}
			}
		} while (avg_diff < -hysis);

		caculate_position(buf, st[0], len[0], pos);
		h_iter->point = pos;
		//printf("point(%d.%d) = %d\n", pos.index, pos.index_decimal, pos.val);
	}
}

void debug_output_peak_vector(vector<struct peak> &vector_peak)
{
	for (vector<struct peak>::iterator iter = vector_peak.begin();
		iter != vector_peak.end(); iter++) {
		printf("%d: %d ", iter->index, iter->val);
		if (iter->sign < 0)
			printf("MAX");
		else if (iter->sign > 0)
			printf("MIN");
		else
			printf("--");
		printf("\n");
	}

	getchar();
}

void debug_output_peak_hamp(vector<struct hamp> &vector_hamp)
{
	for (vector<struct hamp>::iterator iter = vector_hamp.begin();
		iter != vector_hamp.end(); iter++) {
		printf("%d (%d.%d), [%d] %d ~ [%d] %d\n", 
			iter->point.val, iter->point.index, iter->point.index_decimal, 
			iter->full[START]->index, iter->full[START]->val, 
			iter->full[END]->index,iter->full[END]->val);
	}

	getchar();
}

int _tmain(int argc, _TCHAR* argv[])
{
	vector<struct peak> vector_peak;
	vector<struct hamp> vector_hamp;
	int num;
	int a[] = {16,15,-1, 2, 2, 5, 7,10,12,15,13,8,7,-1, 2, 3, -4,  2, 2 , 1, 11, 4, 1, 8, 7, 9,12, 23 ,24, 25, 23,24,24,24,24,24,23,24,23,24,25,20,21};
	int thld = 5;
	int hysis = 2;
	struct merge_config merge_conf = {5, 3, 2};

	num = sizeof(a) / sizeof(a[0]);
	for (int m = 0; m < num; m++)
		cout << a[m] << " ";
	cout << endl;

	findPeaks(a, num, hysis, vector_peak);
	//debug_output_peak_vector(vector_peak);
	
	append_peaks_in_edge(a, num, vector_peak, thld);
	//debug_output_peak_vector(vector_peak);

	merge_adjacent_high_peak(vector_peak, thld, merge_conf, vector_hamp);
	//debug_output_peak_hamp(vector_hamp);

	caculate_peak_in_hamp(a, vector_hamp, thld, hysis);
	debug_output_peak_hamp(vector_hamp);

	return 0;
}
