import scipy
import matplotlib.pyplot as plt
import argparse
from argparse import HelpFormatter



def main():

	parser = argparse.ArgumentParser(prog='VISOR', description='''Script to plot precision, recall and F1 scores''', epilog='''This script has been modified by Davide Bolognin at EMBL/EBI''', formatter_class=CustomFormat) 

	required = parser.add_argument_group('Required I/O arguments')

	required.add_argument('-p','--points', help='Couples of points (precision-recall) to plot', dest='points', type=points, nargs='+', metavar='', required=True)
	required.add_argument('-c', '--colors', help='One color for each couple of precision-recall points', metavar='', nargs='+', action='append', required=True)
	required.add_argument('-m', '--markers', help='One matplotlib marker for each couple of precision-recall points', metavar='', nargs='+', action='append', required=True)
	required.add_argument('-l', '--labels', help='One label for each couple of precision-recall points', metavar='', nargs='+', action='append', required=True)
	required.add_argument('-t', '--title', help='Title of the plot', metavar='', required=True)
	required.add_argument('-lt', '--legendtitle', help='Title of the legend', metavar='', required=True)


	args = parser.parse_args()

	poi=args.points
	colors=args.colors[0]
	markers=args.markers[0]
	labels=args.labels[0]
	title=str(args.title)
	legendtitle=str(args.legendtitle)

	plotPrecisionRecallDiagram(poi, colors, markers, labels, title,legendtitle)

	plt.show()





def points(s):

    try:

        x, y= map(float, s.split(','))

        return x,y

    except:

        raise argparse.ArgumentTypeError('points must be precision,recall')





def fmeasure(p, r):

	"""Calculates the f1 score for precision p and recall r."""

	return float(2*p*r / (p+r))


def fmeasureCurve(f, p):

	"""For a given f1 value and precision get the recall value."""

	return float(f * p / (2 * p - f))


def plotFMeasures(fstepsize=.1, stepsize=0.001):

	"""Plots 10 f1 score curves."""

	p = scipy.arange(0., 1., stepsize)[1:]

	for f in scipy.arange(0., 1., fstepsize)[1:]:

		points = [(x, fmeasureCurve(f, x)) for x in p if 0 < fmeasureCurve(f, x) <= 1.5]
		xs, ys = zip(*points)
		curve, = plt.plot(xs, ys, "--", color="lightgray", linewidth=0.5)
		plt.annotate(r"$f=%.1f$" % f, xy=(xs[-10], ys[-10]), xytext=(xs[-10] - 0.05, ys[-10] - 0.035), size="small", color="black")




def plotPrecisionRecallDiagram(points, colors, markers, labels, title, legendtitle):

	ax = plt.gca()   
	plt.title(title)
	plt.xlabel("Precision")
	plt.ylabel("Recall")
	plotFMeasures()
	scps = []
	plt.axis([-0.02, 1.02, -0.02, 1.02])

	
	for h, (x,y) in enumerate(points):
			
		label = labels[h]
		color=colors[h]
		marker = markers[h]
		scp = ax.scatter(x, y, label=label, c=color,marker=marker, s=30, alpha=0.5)
		scps.append(scp)
				
	plt.legend(loc='lower left', scatterpoints=1, numpoints=1, fancybox=True, title=legendtitle) #loc=0 guess the best location for legend



class CustomFormat(HelpFormatter):

	def _format_action_invocation(self, action):

		if not action.option_strings:

			default = self._get_default_metavar_for_positional(action)
			metavar, = self._metavar_formatter(action, default)(1)
			
			return metavar

		else:

			parts = []

			if action.nargs == 0:

				parts.extend(action.option_strings)

			else:

				default = self._get_default_metavar_for_optional(action)
				args_string = self._format_args(action, default)
				
				for option_string in action.option_strings:

					parts.append(option_string)

				return '%s %s' % (', '.join(parts), args_string)

			return ', '.join(parts)

	def _get_default_metavar_for_optional(self, action):

		return action.dest.upper()






if __name__ == '__main__':

	main()

