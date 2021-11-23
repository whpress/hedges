/* usage:
HeapScheduler<Doub,CargoClass> myheap;
myheap.push(time, cargo);
timeval = myheap.pop(cargoval);  // cargo is returned in cargoval
// or
timeval = myheap.pop();
cargoval = myheap.lastcargo;
// or, if no cargo, can just do
HeapScheduler<> myheap;  // (equivalent to <Doub,void*>, actually)
myheap.push(time);
timeval = myheap.pop();
*/

template <class T = Doub, class U = void *>
struct HeapScheduler
{
	static const Int defaultps = 1100000; // initial heap size
	T bigval;
	U lastcargo;
	Int ps, ks;
	NRvector<T> ar; // times
	NRvector<U> br; // "cargo"

	HeapScheduler() : bigval(numeric_limits<T>::max()), ps(defaultps), ks(0), ar(ps, bigval), br(ps) {}
	void push(T time, U cargo = U(NULL))
	{	// lengthen list, add to end, sift up
		// pushes a time and cargo onto the heap
		Int k, mo;
		if (ks == ps)
			resizear(2 * ps);
		k = ks++;
		ar[k] = time;
		br[k] = cargo;
		while (k > 0 && ar[mo = (k - 1) / 2] > ar[k])
		{
			SWAP(ar[k], ar[mo]); // swap with mother
			SWAP(br[k], br[mo]);
			k = mo;
		}
	}
	T pop() { return pop(lastcargo); } // if no argument, return cargo in HeapScheduler::lastcargo
	T pop(U &cargo)
	{	// return top of heap, move last to top, shorten list, sift down
		// pops the next (in order) time and its cargo from the heap
		// returns numeric_limits::max() time, and U() cargo, when heap is empty
		Int k = 0, rdau, ldau, mindau;
		T ans = ar[0];
		U cans = br[0];
		if ((ks--) > 0)
		{
			ar[0] = ar[ks];
			br[0] = br[ks];
			ar[ks] = bigval;
			br[ks] = U();
			while ((ldau = 2 * k + 1) < ks)
			{
				rdau = ldau + 1; // might be a bigval, but that is OK
				mindau = (ar[ldau] < ar[rdau] ? ldau : rdau);
				if (ar[k] > ar[mindau])
				{
					SWAP(ar[k], ar[mindau]); // swap with smaller of two daughters
					SWAP(br[k], br[mindau]);
				}
				else
					break;
				k = mindau;
			}
		}
		cargo = cans;
		return ans;
	}
	void resizear(Int newps)
	{							// only used internally
		ar.resize(newps, true); //resize preserving contents
		br.resize(newps, true);
		for (int i = ks; i < newps; i++)
			ar[i] = bigval;
		ps = newps;
	}
	void rewind()
	{ // zero out the heap w/o changing its size in memory
		ks = 0;
		ar.assign(ps, bigval);
		br.resize(ps);
		//for (int i = 0; i < ps; i++) {
		//	ar[i] = bigval;
		//	br[i] = U(NULL);
		//}
	}
	void reinit()
	{ // zero out the heap and give back memory
		ks = 0;
		ps = defaultps;
		ar.assign(ps, bigval);
		br.resize(ps);
	}

	/*	void printheap() { // only used for debugging
		printf("heap is:\n");
		for (int i=0;i<ks;i++) printf("%d  %.5f  %.5f\n",i,ar[i],br[i]);
		printf("end heap\n");
	}
*/
};
