/**
 * @fileOverview Custom JQuery code loaded with html page.
 * @see 'sam/browse.html'
 * @author Jeffrey Zhou
 * @version 0.0.1
 */


 /*****************************************************************************
  * Main Method
  ****************************************************************************/
$(function() {

  // Toggle form mode between range/custom
  $('#coordinateForm .btn-group label').click(function(event) {
    var mode = $(this).attr('for');
    if (mode == 'modeInput-range') {
      hideDivsByAttr('data-toggle', 'custom-toggle');
      $('div[data-toggle=range-toggle]').show();
    }
    if (mode == 'modeInput-custom') {
      hideDivsByAttr('data-toggle', 'range-toggle');
      $('div[data-toggle=custom-toggle]').show();
    }
  });

  // Select a Gene from the GeneTable and add it to the form
  $('#selectGeneTable tbody tr').click(function(event) {
    resetGeneInput();
    $(this).addClass('bg-dark text-white');
    $('#geneInput').val($(this).attr('data-gene'));
  });

});


/*****************************************************************************
 * Miscellaneous Functions
 ****************************************************************************/

/**
 * Empty the #geneInput control and clear the #selectGeneTable
 */
var resetGeneInput = function() {
  $('#selectGeneTable tbody tr').removeClass('bg-dark text-white');
  $('#geneInput').val('');
}

/**
 * Hide and reset divs in a form by an attribute-value pair
 * @param {string} attribute The div attribute to select.
 * @param {string} value The div attribute value to match.
 */
var hideDivsByAttr = function(attribute, value) {
  var divSelector = 'div[' + attribute + '=' + value + ']';
  var inputSelector = divSelector + ' input,textarea';
  $(inputSelector).val('');
  $(divSelector).hide();
}
