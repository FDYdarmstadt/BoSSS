preparation:
>>>>>>>>>insert into controlfile>>>>>>>>>>>>>>>
...
//Partitioning
c.GridPartType = GridPartType.Predefined;
c.GridPartOptions = "hallo";

Func<double[], int> MakeMyPartioning = delegate (double[] X) {
	double x = X[0];
	double y = X[1];

	double[] separation = new double[] { 1, 1 };
	switch (cores) {
		case 4:
			separation = new double[] { 2, 2 };
			break;
		case 8:
			separation = new double[] { 4, 2 };
			break;
		case 16:
			separation = new double[] { 4, 4 };
			break;
		case 32:
			separation = new double[] { 8, 4 };
			break;
		case 64:
			separation = new double[] { 8, 8 };
			break;
		default:
			c.GridPartType = GridPartType.none;
			break;
	}

	double xspan = (xMax - xMin) / separation[0];
	double yspan = (yMax - yMin) / separation[1];
	int rank=int.MaxValue;
	int icore = 0;
	for(int i=0; i < separation[0]; i++) {
		for (int j = 0; j < separation[1]; j++) {
			bool xtrue =  x <= xspan * (i + 1)+xMin;
			bool ytrue =  y <= yspan * (j + 1)+yMin;
			if (xtrue && ytrue) {
				rank = icore;
				return rank;
			}
			icore++;
		}
	}

	return rank;
};
...
grid.AddPredefinedPartitioning("hallo", MakeMyPartioning);
...

>>>>>>>>>>>>>>>insert into Main/RunSolverOneStep>>>>>>>>>>>>>>>>
...
using (var ht = new FuncTrace()) {
...
	if (TimestepNo <= 5)
		ilPSP.Tracing.Tracer.Current.ResetRecursive();
}
...

execution:
1.part1
2.part2
3.done