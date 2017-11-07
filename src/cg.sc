import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.extensions.gatk.CombineGVCFs
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.ReferenceConfidenceMode
import org.broadinstitute.gatk.utils.commandline.Output
import org.broadinstitute.gatk.utils.variant.GATKVCFIndexType

class CombineOurGs extends QScript {
  @Input(doc="The input bam file", shortName = "I")
  var bamFile: List[File] = Nil

  @Input(doc="The reference", shortName="R")
  var referenceFile: File = _

  @Output(doc="output", shortName="O")
  var outFile: File = _

  def script(): Unit = {
    val genotyper = new CombineGVCFs
    genotyper.scatterCount = 200
    genotyper.variant = this.bamFile
    genotyper.reference_sequence = this.referenceFile
    genotyper.out = outFile
    add(genotyper)

  }

}