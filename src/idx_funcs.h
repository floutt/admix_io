#include <dotgeno.h>

typedef struct {
	struct idx_head** elems;
	size_t length;
} idx_list_arr;

idx_list_arr init_idx_list_arr(size_t length) {
	idx_list_arr out;
	out.length = length;
	out.elems = (struct idx_head**)malloc(length * sizeof(struct idx_head*));
	return out;
}

void free_idx_list_arr(idx_list_arr* ila) {
	free(ila->elems);
}

void insert_range(size_t length, struct idx_head* head) {
	if(!TAILQ_EMPTY(head)) {
		fprintf(stderr, "ERROR: list must be empty to run insert_range function.\n");
		exit(EXIT_FAILURE);
	}
	for(size_t i = 0; i < length; i++) {
		struct idx_node* idn = (struct idx_node*)malloc(sizeof(struct idx_node));
		idn->idx = i;
		TAILQ_INSERT_TAIL(head, idn, nodes);		
	}
}

/* sorted list processing functions */
size_t get_max_index(struct idx_node** arr, size_t length) {
        size_t max = 0;
        size_t max_idx = 0;
        for(size_t i = 0; i < length; i++) {
                if(arr[i]->idx >= max) {
                        max = arr[i]->idx;
                        max_idx = i;
                }
        }
        return max_idx;
}

bool all_equal(struct idx_node** arr, size_t length) {
        size_t first = arr[0]->idx;
        for(size_t i = 1; i < length; i++) {
                if(first != arr[i]->idx) {
                        return false;
                }
        }
        return true; 
}

bool all_null(struct idx_node** arr, size_t length) {
        for(size_t i = 0; i < length; i++) {
			if(arr[i]) {
				return false;
			}
        }
        return true; 
}
size_t intersect_idx(idx_list_arr* ila, struct idx_head* head_out) {
	size_t cnt = 0;
	// get first elements
	struct idx_node** cur_elems = (struct idx_node**)malloc(ila->length * sizeof(struct idx_node*));
	for(size_t i = 0; i < ila->length; i++) {
		if(TAILQ_EMPTY(ila->elems[i])) { free(cur_elems); return cnt; }
		cur_elems[i] = TAILQ_FIRST(ila->elems[i]);
	}
	
	// fill in output linked list with elements of first list
	struct idx_node* tmp_node;
	TAILQ_FOREACH(tmp_node, ila->elems[0], nodes) {
		struct idx_node* idn = (struct idx_node*)malloc(sizeof(struct idx_node));
		idn->idx = tmp_node->idx;
		TAILQ_INSERT_TAIL(head_out, idn, nodes);
	}
	struct idx_node* cur_elem_out = TAILQ_FIRST(head_out);

	while(true) {
		bool end_while = false;  // change value to exit while loop	
		if(all_equal(cur_elems, ila->length)) {
			cnt++;
			cur_elem_out = TAILQ_NEXT(cur_elem_out, nodes);
			for(size_t i = 0; i < ila->length; i++) {
				cur_elems[i] = TAILQ_NEXT(cur_elems[i], nodes);
				if((cur_elems[i] == NULL)) {
					if(i == 0) { end_while = true; break; }
					struct idx_node* old_val;
					while(cur_elem_out) {
						old_val = cur_elem_out;
						cur_elem_out = TAILQ_NEXT(cur_elem_out, nodes);
						TAILQ_REMOVE(head_out, old_val, nodes);
						free(old_val);
					}
					end_while = true;
					break;
				}
			}
		} else {
			size_t max_i = get_max_index(cur_elems, ila->length);
			for(size_t i = 0; i < ila->length; i++) {
				if(i == max_i) { continue; }
				if(cur_elems[i]->idx == cur_elems[max_i]->idx) { continue; }
				if(i == 0) {
					struct idx_node* old_val = cur_elem_out;
					cur_elem_out = TAILQ_NEXT(cur_elem_out, nodes);
					TAILQ_REMOVE(head_out, old_val, nodes);
					free(old_val);
				}
				cur_elems[i] = TAILQ_NEXT(cur_elems[i], nodes);
				if(cur_elems[i] == NULL) {
					if(i == 0) {
						end_while = true;
						break;
					} else {
						free_idx_list(head_out);
						end_while = true;
						break;
					}
				}
			}
		}
		if(end_while) {
			free(cur_elems);
			break;
		}
	}
	return cnt;
}

void remove_idx(idx_list_arr* ila, struct idx_head* ref_head) {
	struct idx_node** cur_elems_ila = (struct idx_node**)malloc(ila->length * sizeof(struct idx_node*));
	for(size_t i = 0; i < ila->length; i++) {
		if(TAILQ_EMPTY(ila->elems[i])) { return; }
		cur_elems_ila[i] = TAILQ_FIRST(ila->elems[i]);
	}
	struct idx_node* cur_elem_ref = TAILQ_FIRST(ref_head);

	while(1) {
		if(cur_elem_ref == NULL) { break; }
		if(all_null(cur_elems_ila, ila->length)) { break; }
		bool cur_elem_removed = false;
		for(size_t i = 0; i < ila->length; i++) {
			if(cur_elems_ila[i] == NULL) { continue; }
			// if any elements of the remove array are less than the current value in the ref then skip
			// until that is the case
			while(cur_elems_ila[i]->idx < cur_elem_ref->idx) {
				cur_elems_ila[i] = TAILQ_NEXT(cur_elems_ila[i], nodes);
				if(cur_elems_ila[i] == NULL) { break; }
			}
			// if the element is equal to the current value then remove it and iterate
			if(cur_elems_ila[i]) {
				if(cur_elems_ila[i]->idx == cur_elem_ref->idx) {
					struct idx_node* old_val = cur_elem_ref;
					TAILQ_REMOVE(ref_head, cur_elem_ref, nodes);
					cur_elem_ref = TAILQ_NEXT(cur_elem_ref, nodes);
					free(old_val);
					cur_elem_removed = true;
					break;
				}
			}
		}
		// if no element was removed then go to the next element
		if(!cur_elem_removed) { cur_elem_ref = TAILQ_NEXT(cur_elem_ref, nodes); }
	}
	free(cur_elems_ila);
}
